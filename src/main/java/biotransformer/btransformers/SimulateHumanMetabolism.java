package biotransformer.btransformers;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.esaprediction.ESSpecificityPredictor;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.HumanMetabolismBiotransformerFunctions;
import biotransformer.utils.HumanMetabolismHelpingFunctions;
import biotransformer.utils.Utilities;
import biotransformer.validateModels.MultivariateDistribution;
//import reactantpredictor.BioTransfor
import biotransformerapi.BioTransformerAPI;

public class SimulateHumanMetabolism {
	HumanMetabolismHelpingFunctions hm_helper = new HumanMetabolismHelpingFunctions();
	HumanMetabolismBiotransformerFunctions hm_transform = new HumanMetabolismBiotransformerFunctions();
	Double scoreThreshold = 0.0;
	//BioTransformerAPIs reactantPred = new BioTransformerAPIs();
	MultivariateDistribution mnd;
	
	protected LinkedHashMap<String, MetabolicReaction> combinedReactionsHash = new LinkedHashMap<String, MetabolicReaction>();
	
	Cyp450BTransformer cyb = new Cyp450BTransformer(BioSystemName.HUMAN);
	Phase2BTransformer p2b = new Phase2BTransformer(BioSystemName.HUMAN);
	ECBasedBTransformer ecb	= new ECBasedBTransformer(BioSystemName.HUMAN);
	HGutBTransformer hgb = new HGutBTransformer();
	int cyp450Mode = 3;
	
	public  SimulateHumanMetabolism(Integer cyp450Mode) throws Exception{
		this.mnd =  new MultivariateDistribution();
		this.cyp450Mode = cyp450Mode;
		initializeReactionHash();
	}
	/**
	 * This function will simulate the human metabolism process for a small molecule taken orally or intravenously 
	 * The metabolism process is:
	 * GutMicrobial metabolism ==> EC-Commision metabolism ==> Liver Metabolism
	 * 	Here Liver Metabolism is running PhaseI ==> Phase II metabolism 3 iterations
	 * 		Iteration 1: A molecule is goes to Phase I metabolism then itself and its products go to Phase II metabolism, the predicted metabolites are then stored in the metabolites Set
	 * 		Iteration 2: the molecule and its non-Phase II metabolites go to Phase I and themselves with corresponding metabolites go to Phase II again.
	 * 		Iteration 3: repeat the above process one more time
	 * 	Note that a liver metabolite is validated by a multivariate function learned on a dataset of urine metabolites.
	 * 	Please also note that because the small molecule bile metabolites are produced through liver conjugated metbaolites, we don't consider fecal metabolites in the training dataset of the aforementioned multi-variate mdoel
	 *  In addition, the phase II metabolites will be stored into the metabolites set directly and will not have another iteration of CYP450/Phase II metabolism
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	
	public ArrayList<Biotransformation> simulateHumanMetabolism(IAtomContainer oneMole) throws Exception{
		ArrayList<Biotransformation> result = new ArrayList<>();
		ArrayList<Biotransformation> biotransformations_gut = this.hm_transform.runGutMicrobialTransformation_multiple(oneMole, 4);
							
		return result;
	}
	/**
	 * This function will simulate the human liver metabolism only by running a combination of PhaseI and PhaseII metabolism
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> simulateLiverMetabolism(IAtomContainer oneMole) throws Exception{
		ArrayList<Biotransformation> resultBiotransformationList = runPhaseIandPhaseII_byIterations(oneMole, 3);
		return resultBiotransformationList;
	}
	
	public ArrayList<Biotransformation> runPhaseIandPhaseII_byIterations(IAtomContainer startingCompound,  int  numIteration) throws Exception{		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet currentMetabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet startingCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		startingCompounds.addAtomContainer(startingCompound);
		currentMetabolites.addAtomContainer(startingCompound);
		ArrayList<String> subStrateProcessed = new ArrayList<>();
		System.out.println(this.mnd.sg.create(startingCompound));
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		String subStrate_inChiKey = inchiGenFactory.getInChIGenerator(startingCompound).getInchiKey();
		subStrateProcessed.add(subStrate_inChiKey);
		Biotransformer bf = new Biotransformer(BioSystemName.HUMAN);		
		//If a molecule has been used as substrate and it occurs as metablite in another biotransformation, we should not use it as a substrate in the next iteration
		ArrayList<String> processedMetaboliteList = new ArrayList<>();
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		String inchiKey_start = inchiFactory.getInChIGenerator(startingCompound).getInchiKey();
		processedMetaboliteList.add(inchiKey_start);
		for(int i = 0; i < numIteration; i++){
			ArrayList<Biotransformation> result = runPhaseIPhaseIIParallelly(startingCompounds);
			biotransformations.addAll(result);
			ArrayList<Biotransformation> biotransformations_temp = this.mnd.extractValidBioTransformation_Liver(biotransformations);
			IAtomContainerSet nextSubstrates = bf.extractProductsFromBiotransformations(biotransformations_temp);
			IAtomContainerSet nextNewSubstrates = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			for(int j = 0; j < nextSubstrates.getAtomContainerCount(); j++){
				subStrate_inChiKey = inchiGenFactory.getInChIGenerator(nextSubstrates.getAtomContainer(j)).getInchiKey();
				if(!subStrateProcessed.contains(subStrate_inChiKey)){
					subStrateProcessed.add(subStrate_inChiKey);
					try{
						BioTransformerAPI reactantPred = new BioTransformerAPI();
						boolean validCyp450 = reactantPred.predictReactant(nextSubstrates.getAtomContainer(j), ESSpecificityPredictor.cyp450EnzymeList);
						if(!validCyp450) continue;
					}catch (Exception e){
						nextNewSubstrates.addAtomContainer(nextSubstrates.getAtomContainer(j));
					}
					nextNewSubstrates.addAtomContainer(nextSubstrates.getAtomContainer(j));
				}
			}
			//System.out.println("--------------" + nextNewSubstrates.getAtomContainerCount() + "---------------");			
			startingCompounds = getUniqueMetabolites(nextNewSubstrates, processedMetaboliteList);
					
		}
		biotransformations = Utilities.selectUniqueMetabolites(biotransformations);
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	/**
	 * Special Sequence: Run Phase I and Phase II parallelly
	 * This function will run one Phase I metabolism prediction and Phse II metabolism prediction and merge the metabolites.
	 * The merged metabolites are returned as metabolites.
	 */
	public ArrayList<Biotransformation> runPhaseIPhaseIIParallelly(IAtomContainerSet startingCompounds) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		LinkedHashMap<Biotransformer.bType, Object> btransformers = new LinkedHashMap<Biotransformer.bType, Object>();
		IAtomContainerSet currentMetabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<Biotransformation> currentBiots =  new ArrayList<Biotransformation>();
		//currentMetabolites.addAtomContainer(startingCompound);
		/**
		 * CYP450
		 */
	
		IAtomContainerSet cyp450Substrates = (IAtomContainerSet) startingCompounds.clone();
		//System.out.println("-------------- CYP450 substrates" + cyp450Substrates.getAtomContainerCount() + "---------------");
		//btransformers.put(bType.CYP450, new Cyp450BTransformer(BioSystemName.HUMAN));		
		//currentBiots = ((Cyp450BTransformer) btransformers.get(bType.CYP450)).predictCyp450BiotransformationChain(cyp450Substrates, true, true, 1, this.scoreThreshold);	
		//currentMetabolites.add(((Biotransformer) btransformers.get(bType.CYP450)).extractProductsFromBiotransformations(currentBiots));
		//currentBiots = this.cyb.predictCyp450BiotransformationChain(cyp450Substrates, true, true, 1, this.scoreThreshold);
		currentBiots = this.cyb.predictCyp450BiotransformationsByMode(cyp450Substrates,this.cyp450Mode, true, true, this.scoreThreshold);
		currentMetabolites.add(this.cyb.extractProductsFromBiotransformations(currentBiots));
		biotransformations.addAll(currentBiots);
		/**
		 * Phase II
		 */
		IAtomContainerSet phaseIISubstrates = (IAtomContainerSet) startingCompounds.clone();
		//System.out.println("-------------- PhseII substrates" + phaseIISubstrates.getAtomContainerCount() + "---------------");
		//btransformers.put(bType.PHASEII, new Phase2BTransformer(BioSystemName.HUMAN));		
		//currentBiots = ((Phase2BTransformer) btransformers.get(bType.PHASEII)).applyPhase2TransformationsChainAndReturnBiotransformations(phaseIISubstrates, true, true, true, 1, this.scoreThreshold);		
		//currentMetabolites.add(((Biotransformer) btransformers.get(bType.PHASEII)).extractProductsFromBiotransformations(currentBiots));
		currentBiots = this.p2b.applyPhase2TransformationsChainAndReturnBiotransformations(phaseIISubstrates, true, true, true, 1, this.scoreThreshold);
		currentMetabolites.add(this.p2b.extractProductsFromBiotransformations(currentBiots));
		biotransformations.addAll(currentBiots);
		biotransformations = Utilities.selectUniqueMetabolites(biotransformations);
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	/**
	 * This function will extract unique metabolites from all molecules in the query molecules stored in the IAtomContainerSet
	 * @param allMolecules
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet getUniqueMetabolites(IAtomContainerSet allMolecules, ArrayList<String> processedList) throws Exception{
		IAtomContainerSet cleanedSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		for(int i = 0; i < allMolecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = allMolecules.getAtomContainer(i);
			if(isValidMetabolte(oneMole)) continue;
			InChIGenerator inchiGen = inchiFactory.getInChIGenerator(oneMole);
			String inChiKey = inchiGen.getInchiKey();
			if(!processedList.contains(inChiKey)){
				processedList.add(inChiKey);
				cleanedSet.addAtomContainer(oneMole);
			}
		}
		return cleanedSet;
	}
	/**
	 * This function will initialize the combinedReactionsHash by adding all reactions from all biotransformer types
	 */
	public void initializeReactionHash(){
		for(Map.Entry<String, MetabolicReaction> m : this.ecb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(m.getKey())){
				this.combinedReactionsHash.put(m.getKey(), m.getValue());
			}
		}		
		for(Map.Entry<String, MetabolicReaction> n : this.cyb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(n.getKey())){
				this.combinedReactionsHash.put(n.getKey(), n.getValue());
			}
		}
		
		for(Map.Entry<String, MetabolicReaction> p : this.hgb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(p.getKey())){
				this.combinedReactionsHash.put(p.getKey(), p.getValue());
			}
		}		
		for(Map.Entry<String, MetabolicReaction> o : this.p2b.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(o.getKey())){
				this.combinedReactionsHash.put(o.getKey(), o.getValue());
			}
		}
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
	/**
	 * This function will save the predicted results into sdf file
	 * @param molecule
	 * @param outputFileName
	 * @param annotate
	 * @throws Exception
	 */
	public void simulateHumanMetabolismAndSaveToSDF(IAtomContainer molecule, String outputFileName, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations = this.simulateLiverMetabolism(molecule);
		this.cyb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}

}
