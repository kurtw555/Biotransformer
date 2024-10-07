package biotransformer.validateModels;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.Utilities;
//import predictor.PredictLogD;

public class MultivariateDistribution {
	public SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	public SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
	public double threadhold_prob;
	public MultivariateNormalDistribution mnd;
	public HashMap<String, HashMap<Double,Double>> featureValueProbMap;
	public FeatureGeneration mvd = new FeatureGeneration();
	public String[] titleLineList;
//	public PredictLogD logD_predictor = new PredictLogD();
	
	public static void main(String[] args) throws Exception{
		MultivariateDistribution mnd = new MultivariateDistribution("GutMetabolite");
		String smiles = "OC1=CC(O)=C2C(=O)C=C(OC2=C1)C1=CC(O)=C(O)C=C1";
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer oneMole = sp.parseSmiles(smiles);
		System.out.println(mnd.makePrediction(oneMole, "GutMetabolite"));
	}

	
	public MultivariateDistribution(String type) throws Exception{
		String mnd_path;
		String discreateFeatuerPath;
		if(type.equalsIgnoreCase("Liver")){
			mnd_path = "./supportfiles/MVDModels/LiverModelParameters/mnd_model_parameters.csv";		
			//discreateFeatuerPath = "./supportfiles/MVDModels/LiverModelParameters/bayesian_model_parameters.csv";
			discreateFeatuerPath = null;
			this.threadhold_prob = 1.5E-43;
		}
		else if(type.equalsIgnoreCase("ECBased")){
			mnd_path = "./supportfiles/MVDModels/ECBasedModelParameters/mnd_model_parameters.csv";		
			//discreateFeatuerPath = "./supportfiles/MVDModels/ECBasedModelParameters/bayesian_model_parameters.csv";
			discreateFeatuerPath = null;
			this.threadhold_prob = 1.5E-43;
		}
		else if(type.equalsIgnoreCase("GutSubstrate")){
			mnd_path = "./supportfiles/MVDModels/GutSubstratesModelParameters/mnd_model_parameters.csv";		
			discreateFeatuerPath = "./supportfiles/MVDModels/GutSubstratesModelParameters/bayesian_model_parameters.csv";
			this.threadhold_prob = 1.0E-39;
		}
		else if(type.equalsIgnoreCase("GutMetabolite")){
			mnd_path = "./supportfiles/MVDModels/GutMetabolitesModelParameters/mnd_model_parameters.csv";		
			discreateFeatuerPath = "./supportfiles/MVDModels/GutMetabolitesModelParameters/bayesian_model_parameters.csv";
			this.threadhold_prob = 1.0E-43;
		}
		else throw new Exception("The multivariateDistruction can not be initialzed on a unknown Metabolism type which is " + type.toUpperCase());
		this.featureValueProbMap = readDataForBayes(discreateFeatuerPath);
		getMultivariateDistModel(mnd_path);
		this.titleLineList = this.mvd.generateTitleLine().split(",");
	}
	
	/**
	 * This function will predict the (probability density of continuous features * probability of functional discrete features) of a 
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public double makePrediction(IAtomContainer oneMole, String type) throws Exception{
		String smiles = this.sg.create(oneMole);
		//System.out.println(smiles);
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
		String[] targetFeaturesToRemove = null;
		if(type.equals("Liver")) targetFeaturesToRemove = FeatureGeneration.featuresToRemove_Urine;
		else if(type.equals("ECBased")) targetFeaturesToRemove = FeatureGeneration.featuresToRemove_ECBased;
		else if(type.equals("GutMetabolite")) targetFeaturesToRemove = FeatureGeneration.featuresToRemove_gut_metabolite;
		else if(type.equals("GutSubstrate")) targetFeaturesToRemove = FeatureGeneration.featuresToRemove_gut_substrate;
		else throw new Exception("Can't figure out the metbaolism type");
		for(int i = 0; i < this.titleLineList.length; i++){
			boolean isValid = true;
			String oneRawFeature = this.titleLineList[i];		
			for(int j = 0; j < targetFeaturesToRemove.length; j++){
				if(oneRawFeature.equals(targetFeaturesToRemove[j])){
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
		if(this.featureValueProbMap == null || this.featureValueProbMap.isEmpty()) return resultProb;
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
		if(inputPath == null || inputPath.equals("")) return null;
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
		br.close();
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
	public ArrayList<Biotransformation> extractValidBioTransformation_CYP450(ArrayList<Biotransformation> transformationList) throws Exception{		
		
		ArrayList<Biotransformation> resultArray = new ArrayList<>();
		String[] cyp450List = {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4s"};
		
		for(int i = 0; i < transformationList.size(); i++){
			Biotransformation bt = transformationList.get(i);
			String precursor_SMILES = this.sg.create(bt.getSubstrates().getAtomContainer(0));
			IAtomContainer precursor = this.sp.parseSmiles(precursor_SMILES);
			Double proDensity_precursor = makePrediction(precursor, "Liver");			
			//Double logD_precursor = this.logD_predictor.predictLogD(precursor);
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
			boolean havingEnoughMetabolites = false;
			IAtomContainerSet metabolites = bt.getProducts();
			int oxygen_precursor = countOxygen(precursor);
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				IAtomContainer oneMetabolite = metabolites.getAtomContainer(j);
				if(!isValidMetabolte(oneMetabolite)) continue;
				else havingEnoughMetabolites = true;
				double density_metabolite = makePrediction(oneMetabolite, "Liver");	
				int oxygen_metabolite = countOxygen(oneMetabolite);
				//Double logD_metabolite = this.logD_predictor.predictLogD(metabolites.getAtomContainer(j));
				//Removed the condition: logD_metabolite > logD_precursor ||  Not sure if this is true or not
				if(density_metabolite < proDensity_precursor*0.9 ||  density_metabolite < this.threadhold_prob){
					if(oxygen_precursor >= oxygen_metabolite){
						validMetabolite = false;
						break;
					}
				}
			}
			if(!validMetabolite || !havingEnoughMetabolites) continue;
			else resultArray.add(bt);
		}
		resultArray = Utilities.selectUniqueBiotransformations(resultArray);
		return resultArray;
	}
	/**
	 * This function will count the number of oxygen atoms in the molecule
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public int countOxygen(IAtomContainer molecule) throws Exception{
		int count = 0;
		for(int i = 0; i < molecule.getAtomCount(); i++){
			if(molecule.getAtom(i).getSymbol().equalsIgnoreCase("O")) count++;
		}
		return count;
	}
	/**
	 * This function will remove the invalid biotransformations for ECBased metabolism
	 * A biotransformation is invalid if :
	 * 1. the probability density of its metbaolite is lower than the threshold
	 * 2. The probability density of the metabolite is lower than the precursor's 
	 * @param transformationList
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> extractValidBioTransformation_ECBased(ArrayList<Biotransformation> transformationList) throws Exception{
		ArrayList<Biotransformation> resultArray = new ArrayList<>();
		for(int i = 0; i < transformationList.size(); i++){
			Biotransformation bt = transformationList.get(i);
			String precursor_SMILES = this.sg.create(bt.getSubstrates().getAtomContainer(0));
			//IAtomContainer precursor = this.sp.parseSmiles(precursor_SMILES);
			//Double proDensity_precursor = makePrediction(precursor, "ECBased");	
			IAtomContainerSet metabolites = bt.getProducts();
			boolean validMetabolite = true;
			boolean havingEnoughMetabolites = false;
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				IAtomContainer oneMetabolite = metabolites.getAtomContainer(j);
				if(!isValidMetabolte(oneMetabolite)) continue;
				else havingEnoughMetabolites = true;
//				double density_metabolite = makePrediction(oneMetabolite, "ECBased");				
//				//if(density_metabolite < proDensity_precursor || density_metabolite < this.threadhold_prob){
//				if(density_metabolite < this.threadhold_prob){
//					validMetabolite = false;
//					break;
//				}
			}
			if(!validMetabolite || !havingEnoughMetabolites) continue;
			else resultArray.add(bt);
		}
		resultArray = Utilities.selectUniqueBiotransformations(resultArray);
		return resultArray;
	}
	
	/**
	 * This function will remove the invalid biotransformations for Gut Microbial metabolism
	 * A biotransformation is invalid if :
	 * 1. the probability density of its metbaolite is lower than the threshold
	 * 2. The probability density of the metabolite is lower than the precursor's 
	 * @param transformationList
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> extractValidBioTransformation_GutMetabolites(ArrayList<Biotransformation> transformationList) throws Exception{
		ArrayList<Biotransformation> resultArray = new ArrayList<>();
		for(int i = 0; i < transformationList.size(); i++){			
				Biotransformation bt = transformationList.get(i);
				String precursor_SMILES = this.sg.create(bt.getSubstrates().getAtomContainer(0));
				//IAtomContainer precursor = this.sp.parseSmiles(precursor_SMILES);
			try{
				//Double proDensity_precursor = makePrediction(precursor, "GutMetabolite");	
				IAtomContainerSet metabolites = bt.getProducts();
				boolean validMetabolite = true;
				boolean havingEnoughMetabolites = false;
				for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
					IAtomContainer oneMetabolite = metabolites.getAtomContainer(j);
					if(!isValidMetabolte(oneMetabolite)) continue;
					else havingEnoughMetabolites = true;
					double density_metabolite = makePrediction(oneMetabolite, "GutMetabolite");				
					//if(density_metabolite < proDensity_precursor || density_metabolite < this.threadhold_prob){
					if(density_metabolite < this.threadhold_prob){
						validMetabolite = false;
						break;
					}
				}
				if(!validMetabolite || !havingEnoughMetabolites) continue;
				else resultArray.add(bt);
			}catch (Exception e){
				System.out.println("The precursor can not be identified : " + precursor_SMILES);
			}
		}
		resultArray = Utilities.selectUniqueBiotransformations(resultArray);
		return resultArray;
	}
	/**
	 * This function will Simply check if a molecule is a valid gut microbial substrate by
	 * checking if its probability density is lower than the threshold.
	 * Invalid if it's lower
	 * Valid otherwise
	 * @param transformationList
	 * @return
	 * @throws Exception
	 */
	public boolean isValidGutSubstrate(IAtomContainer oneMole) throws Exception{
		Double proDensity_precursor = makePrediction(oneMole.clone(), "GutSubstrate");
		if(proDensity_precursor >= this.threadhold_prob) return true;
		else return false;		
	}
	/**
	 * This function will check the candidate substrates for gut microbial metabolism and extract those valid ones
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet isValidGutSubstrates(IAtomContainerSet molecules) throws Exception{
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			if(isValidGutSubstrate(molecules.getAtomContainer(i))) resultSet.addAtomContainer(molecules.getAtomContainer(i));
		}
		return resultSet;
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
		if(countAtom_nonH >=2) return true;
		else return false;
	}
}
