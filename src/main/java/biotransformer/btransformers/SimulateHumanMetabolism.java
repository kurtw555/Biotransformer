package biotransformer.btransformers;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.HumanMetabolismHelpingFunctions;
import biotransformer.utils.Utilities;
//import biotransformerapicypreact.BioTransformerAPI_cypreact;
import biotransformerapi.BioTransformerAPI;


public class SimulateHumanMetabolism {
	HumanMetabolismHelpingFunctions hm_helper;
	boolean useDB;
	boolean useSubstitution;
	Double scoreThreshold = 0.0;
	//BioTransformerAPI_cypreact reactantPred = new BioTransformerAPI_cypreact();	
	
	
	protected LinkedHashMap<String, MetabolicReaction> combinedReactionsHash = new LinkedHashMap<String, MetabolicReaction>();

	int cyp450Mode = 1;
	
	public  SimulateHumanMetabolism(Integer cyp450Mode, boolean useDB, String source, boolean checkOverflow, int overflow_threshold, boolean useSubstitution) throws Exception{		
		this.useDB = useDB;
		this.useSubstitution = useSubstitution;
		this.cyp450Mode = cyp450Mode;
		this.hm_helper = new HumanMetabolismHelpingFunctions(0.0, useDB, this.cyp450Mode, source, checkOverflow, overflow_threshold, useSubstitution);
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
	
	public ArrayList<Biotransformation> simulateHumanMetabolism(IAtomContainer oneMole, int iteration) throws Exception{
		//System.out.println(this.hm_helper.sg.create(oneMole));
		ArrayList<Biotransformation> result = new ArrayList<>();
		//SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);

		ArrayList<Biotransformation> biotransformations_gut = this.hm_helper.runCompleteIterations(oneMole,iteration);
		result = biotransformations_gut;			
		return result;
	}
	

	/**
	 * This function will initialize the combinedReactionsHash by adding all reactions from all biotransformer types
	 */
	public void initializeReactionHash(){
		for(Map.Entry<String, MetabolicReaction> m : this.hm_helper.ecb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(m.getKey())){
				this.combinedReactionsHash.put(m.getKey(), m.getValue());
			}
		}		
		for(Map.Entry<String, MetabolicReaction> n : this.hm_helper.cyb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(n.getKey())){
				this.combinedReactionsHash.put(n.getKey(), n.getValue());
			}
		}
		
		for(Map.Entry<String, MetabolicReaction> p : this.hm_helper.hgb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(p.getKey())){
				this.combinedReactionsHash.put(p.getKey(), p.getValue());
			}
		}		
		for(Map.Entry<String, MetabolicReaction> o : this.hm_helper.p2b.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(o.getKey())){
				this.combinedReactionsHash.put(o.getKey(), o.getValue());
			}
		}
	}

	/**
	 * This function will save the predicted results into sdf file
	 * @param molecule
	 * @param outputFileName
	 * @param annotate
	 * @throws Exception
	 */
	public void simulateHumanMetabolismAndSaveToSDF(IAtomContainer molecule, String outputFileName, int iteration, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations = this.simulateHumanMetabolism(molecule, iteration);
		this.hm_helper.cyb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}
	
	/**
	 * These functions are used to save the predicted biotransformations into sdf file or csv file.
	 * @param biotransformations
	 * @param outputFileName
	 * @param annotate
	 * @throws Exception
	 */
	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception {
		try{
			this.hm_helper.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}
	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception {
		try{
			this.hm_helper.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}	
	public InChIGeneratorFactory getInChIGenFactory(){
		return this.hm_helper.ecb.inchiGenFactory;
	}
}
