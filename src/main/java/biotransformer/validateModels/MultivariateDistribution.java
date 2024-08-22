package biotransformer.validateModels;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import biotransformer.transformation.Biotransformation;
import biotransformer.utils.Utilities;
import predictor.PredictLogD;

public class MultivariateDistribution {
	public SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	public SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
	public double threadhold_prob = 1.5E-43;
	public MultivariateNormalDistribution mnd;
	public HashMap<String, HashMap<Double,Double>> featureValueProbMap;
	public FeatureGeneration mvd = new FeatureGeneration();
	public String[] titleLineList;
	public PredictLogD logD_predictor = new PredictLogD();
	
	public static void main(String[] args) throws Exception{
		MultivariateDistribution mnd = new MultivariateDistribution();
		String smiles = "NC(N)=NC1=NC(CSCCC(N)=NS(N)(=O)=O)=CS1";
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer oneMole = sp.parseSmiles(smiles);
		System.out.println(mnd.makePrediction(oneMole));
	}

	
	public MultivariateDistribution() throws Exception{
		String mnd_path = "./supportfiles/LiverModelParameters/mnd_model_parameters.csv";
		getMultivariateDistModel(mnd_path);
		String discreateFeatuerPath = "./supportfiles/LiverModelParameters/bayesian_model_parameters.csv";
		this.featureValueProbMap = readDataForBayes(discreateFeatuerPath);
		this.titleLineList = this.mvd.generateTitleLine().split(",");
	}
	
	/**
	 * This function will predict the (probability density of continuous features * probability of functional discrete features) of a 
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public double makePrediction(IAtomContainer oneMole) throws Exception{
		String smiles = this.sg.create(oneMole);
		//Note that the smiles is used as an identifier in the feature line. It will be automatically ignored when generate "useful" features later because "SMILES" is in the featureNotUsedList
		String allFeatures = smiles + "," + this.mvd.generateMolecularFeature(oneMole);
		String[] allFeaturesLIst = allFeatures.split(",");
		//Generate continuous features for the molecule
		double[] continuousFeatureList = new double[FeatureGeneration.continuousFeatures.length];
		for(int i = 0; i < FeatureGeneration.continuousFeatures.length; i++){
			for(int j = 0; j < this.titleLineList.length; j++){
				if(FeatureGeneration.continuousFeatures[i].equals(this.titleLineList[j])){
					continuousFeatureList[i] = Double.parseDouble(allFeaturesLIst[j]);
					break;
				}
			}
		}
		//Generate discrete features for the molecule
		HashMap<String, Double> discreteFeaturesMap = new HashMap<>();				
		for(int i = 0; i < this.titleLineList.length; i++){
			boolean isValid = true;
			String oneRawFeature = this.titleLineList[i];		
			for(int j = 0; j < FeatureGeneration.featuresToRemove.length; j++){
				if(oneRawFeature.equals(FeatureGeneration.featuresToRemove[j])){
					isValid = false;
					break;
				}
			}
			for(int j = 0; j < FeatureGeneration.continuousFeatures.length; j++){
				if(oneRawFeature.equals(FeatureGeneration.continuousFeatures[j])){
					isValid = false;
					break;
				}
			}
			if(isValid){
				discreteFeaturesMap.put(oneRawFeature, Double.parseDouble(allFeaturesLIst[i]));
			}
		}
		Double probDensity = this.mnd.density(continuousFeatureList);
		Double discreteProb = predictDiscreteProb(discreteFeaturesMap);
		return probDensity * discreteProb;
	}
	/**
	 * This function will predict probability for the discrete values for the query compound
	 * @param discreteFeaturesMap
	 * @return
	 */
	public Double predictDiscreteProb(HashMap<String, Double> discreteFeaturesMap){
		Double resultProb = 1.0;
		for(String feature : discreteFeaturesMap.keySet()){
			Double value = discreteFeaturesMap.get(feature);
			Double prob = this.featureValueProbMap.get(feature).get(value);
			//If that discrete value doesn't occur in the training dataset before, then we assume it's invalid to have that value;
			if(prob == null || prob == 0.0){
				prob = 0.0; 
			}
			else resultProb = resultProb * prob;
		}
		return resultProb;
	}
	/**
	 * This function will prepare the multivariate normal distribution model
	 * @param inputPath_mnd
	 * @throws Exception
	 */
	public void getMultivariateDistModel(String inputPath_mnd) throws Exception{		
		BufferedReader br = new BufferedReader(new FileReader(new File(inputPath_mnd)));
		String oneLine = br.readLine();//mean value line
		String[] parseLine = oneLine.split(",");
		double[] meanList = new double[parseLine.length];
		double[][] coVarianceMatrix = new double[parseLine.length][parseLine.length];
		//set up the mean value list
		for(int i = 0; i < meanList.length; i++){
			meanList[i] = Double.parseDouble(parseLine[i]);
		}
		//set up the covariance matrix
		int rowNumber = 0;
		while((oneLine = br.readLine()) != null){
			parseLine = oneLine.split(",");
			for(int i = 0; i < parseLine.length; i++){
				coVarianceMatrix[rowNumber][i] = Double.parseDouble(parseLine[i]);
			}
			rowNumber++;
			
		}
		br.close();

		this.mnd = trainMultivariateDistModel(meanList, coVarianceMatrix);		
	}
	
	public MultivariateNormalDistribution trainMultivariateDistModel(double[] means, double[][] covariances) throws Exception{
		MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(null,means, covariances);
		return mnd;
	}
	/**
	 * This function will read the data with discrete values from the inputPath for the bayesian model
	 * The results are stored in HashMap<Feature, HashMap<DiscreteValue,Probability>>
	 * @param inputPath
	 * @return
	 * @throws Exception
	 */
	public HashMap<String, HashMap<Double,Double>> readDataForBayes(String inputPath) throws Exception{
		HashMap<String, HashMap<Double,Double>> resultMap = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(inputPath));
		String oneLine;
		while((oneLine = br.readLine()) != null){
			String[] parseLine = oneLine.split(",");
			String feature = parseLine[0];
			HashMap<Double,Double> value_Prob_Map = new HashMap<>();
			for(int i = 1; i < parseLine.length; i++){
				String oneEntry = parseLine[i];
				Double value = Double.parseDouble(oneEntry.split(";")[0]);
				Double prob = Double.parseDouble(oneEntry.split(";")[1]);
				value_Prob_Map.put(value, prob);
			}
			resultMap.put(feature, value_Prob_Map);
		}		
		return resultMap;
	}
	/**
	 * This function will remove the invliad biotransformations from the predicted Biotransformations.
	 * A predicted Biotranformation  is valid if the products within the biotransformation are valid 
	 * 1. The metabolites  should be more likely to be a urine/gut metabolite
	 * 2. Here the liver metabolism should make the metabolite more water soluble than the precursor
	 * @param transformationList
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> extractValidBioTransformation_Liver(ArrayList<Biotransformation> transformationList) throws Exception{		
		
		ArrayList<Biotransformation> resultArray = new ArrayList<>();
		String[] cyp450List = {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4s"};
		
		for(int i = 0; i < transformationList.size(); i++){
			Biotransformation bt = transformationList.get(i);
			String precursor_SMILES = this.sg.create(bt.getSubstrates().getAtomContainer(0));
			IAtomContainer precursor = this.sp.parseSmiles(precursor_SMILES);
			Double proDensity_precursor = makePrediction(precursor);			
			Double logD_precursor = this.logD_predictor.predictLogD(precursor);
			ArrayList<String> enzymeArray = bt.getEnzymeNames();
			boolean isCYP450_metabolite = false;
			for(int t = 0; t < enzymeArray.size(); t++){
				if(isCYP450_metabolite) break;
				String enzyme = enzymeArray.get(t);
				for(int k = 0; k < cyp450List.length; k++){
					if(enzyme.contains(cyp450List[k])){
						isCYP450_metabolite = true;
						break;
					}
				}
			}
			if(!isCYP450_metabolite){
				continue;
			}
			boolean validMetabolite = true;
			IAtomContainerSet metabolites = bt.getProducts();
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				double density_metabolite = makePrediction(metabolites.getAtomContainer(j));				
				Double logD_metabolite = this.logD_predictor.predictLogD(metabolites.getAtomContainer(j));
				if(density_metabolite < proDensity_precursor || logD_precursor < logD_metabolite || density_metabolite < this.threadhold_prob){
					validMetabolite = false;
					break;
				}
			}
			boolean havingEnoughMetabolites = false;
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				if(isValidMetabolte(metabolites.getAtomContainer(j))) havingEnoughMetabolites = true;
			}
			if(!validMetabolite || !havingEnoughMetabolites) continue;
			else resultArray.add(bt);
		}
		resultArray = Utilities.selectUniqueBiotransformations(resultArray);
		return resultArray;
	}
	
	/**
	 * Here we don't want to consider the metabolites that have too few atoms in the structure.
	 * @param oneMole
	 * @return
	 */
	public boolean isValidMetabolte(IAtomContainer oneMole){
		int countAtom_nonH = 0;
		for(int t = 0; t < oneMole.getAtomCount(); t++){
			if(!oneMole.getAtom(t).getSymbol().equals("H")) countAtom_nonH++;
		}
		if(countAtom_nonH >=5) return true;
		else return false;
	}
}
