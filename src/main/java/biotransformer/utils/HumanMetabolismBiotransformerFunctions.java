package biotransformer.utils;

import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemicalClassFinder.ChemicalClassName;

public class HumanMetabolismBiotransformerFunctions {
	HumanMetabolismHelpingFunctions hm_helper = new HumanMetabolismHelpingFunctions();
	
	
	/**
	 * This function will run phaseII transformations on all query molecules
	 * @param molecules
	 * @param cyp450Mode
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runPhaseIITransformations(IAtomContainerSet molecules) throws Exception{
		ArrayList<Biotransformation> result = new ArrayList<>();
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = molecules.getAtomContainer(i);
			ArrayList<Biotransformation> result_one = runPhaseIITransformations(oneMole);
			for(int j = 0; j < result_one.size(); j++){
				if(!result.contains(result_one.get(j))){
					result.add(result_one.get(j));
				}
			}
		}
		return result;
	}
	/**
	 * This function will run PhaseII transformations on the query molecule
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runPhaseIITransformations(IAtomContainer oneMole) throws Exception{
		Phase2BTransformer p2 = new Phase2BTransformer(BioSystemName.HUMAN);
		return p2.applyPhase2TransformationsChainAndReturnBiotransformations(oneMole, true, true, true, 1, 0.0);
	}
	/**
	 * This function will run cyp450 transformations on all query molecules
	 * @param molecules
	 * @param cyp450Mode
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runCYP450Transformation(IAtomContainerSet molecules, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> result = new ArrayList<>();
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = molecules.getAtomContainer(i);
			ArrayList<Biotransformation> result_one = runCYP450Transformation(oneMole, cyp450Mode);
			for(int j = 0; j < result_one.size(); j++){
				if(!result.contains(result_one.get(j))){
					result.add(result_one.get(j));
				}
			}
		}
		return result;
	}
	/**
	 * This function will run cyp450 transformation one the query molecule
	 * @param oneMole
	 * @param cyp450Mode
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runCYP450Transformation(IAtomContainer oneMole, int cyp450Mode) throws Exception{
		Cyp450BTransformer cyp450 = new Cyp450BTransformer(BioSystemName.HUMAN);
		return cyp450.predictCyp450BiotransformationsByMode(oneMole, cyp450Mode, true, true, 0.0);
	}
	/**
	 * This function will run GutMicrobialTransformation multiple times depending on the input iteration value 
	 * Specifically, if iteration = -1, it will keep running until the program thinks there are no more novel gut microbial metaboltes
	 * @param oneMole
	 * @param iteration
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runGutMicrobialTransformation_multiple(IAtomContainer oneMole, int iteration) throws Exception{	
		ArrayList<Biotransformation> resultBioTransformationList = new ArrayList<>();
		ArrayList<String> processedSubstrates_SMILES_list = new ArrayList<>();
		IAtomContainerSet nextSubstrates = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		nextSubstrates.addAtomContainer(oneMole);
		for(int i = 0; i < iteration; i++){		
			IAtomContainerSet nextSubstrates_temp = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			ArrayList<Biotransformation> result_oneIteration = runGutMicrobialTransformations(nextSubstrates);
			ArrayList<Biotransformation> result_oneIteration_valid = this.hm_helper.extractValidGutBiotransformations(result_oneIteration);
			for(int t = 0; t < result_oneIteration_valid.size(); t++){
				Biotransformation biotransformation_one = result_oneIteration_valid.get(t);
				IAtomContainerSet nextSubstrates_raw = this.hm_helper.extractProductsFromBioTransformation_one(biotransformation_one);
				for(int j = 0; j < nextSubstrates_raw.getAtomContainerCount(); j++){
					IAtomContainer nextSubstrate_one = nextSubstrates_raw.getAtomContainer(j);
					String smiles_one = this.hm_helper.sg.create(nextSubstrate_one);
					if(!processedSubstrates_SMILES_list.contains(smiles_one)){
						processedSubstrates_SMILES_list.add(smiles_one);
						nextSubstrates_temp.addAtomContainer(nextSubstrate_one);
						resultBioTransformationList.add(biotransformation_one);
					}
				}				
			}
			nextSubstrates = nextSubstrates_temp;
		}
		if(iteration == -1){
			while(nextSubstrates!=null && !nextSubstrates.isEmpty()){
				IAtomContainerSet nextSubstrates_temp = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
				ArrayList<Biotransformation> result_oneIteration = runGutMicrobialTransformations(nextSubstrates);
				ArrayList<Biotransformation> result_oneIteration_valid = this.hm_helper.extractValidGutBiotransformations(result_oneIteration);
				for(int t = 0; t < result_oneIteration_valid.size(); t++){
					Biotransformation biotransformation_one = result_oneIteration_valid.get(t);
					IAtomContainerSet nextSubstrates_raw = this.hm_helper.extractProductsFromBioTransformation_one(biotransformation_one);
					for(int j = 0; j < nextSubstrates_raw.getAtomContainerCount(); j++){
						IAtomContainer nextSubstrate_one = nextSubstrates_raw.getAtomContainer(j);
						String smiles_one = this.hm_helper.sg.create(nextSubstrate_one);
						if(!processedSubstrates_SMILES_list.contains(smiles_one)){
							processedSubstrates_SMILES_list.add(smiles_one);
							nextSubstrates_temp.addAtomContainer(nextSubstrate_one);
							resultBioTransformationList.add(biotransformation_one);
						}
					}				
				}
				nextSubstrates = nextSubstrates_temp;
			}
		}
		return resultBioTransformationList;
	}
	
	/**
	 * This function will run one gutMicrobial transformation on each of the queried molecules
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runGutMicrobialTransformations(IAtomContainerSet molecules) throws Exception{
		ArrayList<Biotransformation> result = new ArrayList<>();
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = molecules.getAtomContainer(i);
			ArrayList<Biotransformation> result_one = runGutMicrobialTransformations(oneMole);
			for(int j = 0; j < result_one.size(); j++){
				if(!result.contains(result_one.get(j))){
					result.add(result_one.get(j));
				}
			}
		}
		return result;
	}
	/**
	 * This function will run one gutMicrobial transformation on the query onemole
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runGutMicrobialTransformations(IAtomContainer oneMole) throws Exception{
		HGutBTransformer hgut = new HGutBTransformer();		
		//ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(oneMole);
		ArrayList<Biotransformation> biotransformations = hgut.simulateGutMicrobialMetabolism(oneMole, true, true, 1, 0.0);
		return biotransformations;		
	}
	/**
	 * This function will run ECBased transformation on each of the query molecules
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runECBasedTransformation(IAtomContainerSet molecules) throws Exception{
		ArrayList<Biotransformation> result = new ArrayList<>();
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = molecules.getAtomContainer(i);
			ArrayList<Biotransformation> result_one = runECBasedTransformations(oneMole);
			for(int j = 0; j < result_one.size(); j++){
				if(!result.contains(result_one.get(j))){
					result.add(result_one.get(j));
				}
			}
		}
		return result;
	}
	/**
	 * This function will run one ECBasedTransformation on one query molecule
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runECBasedTransformations(IAtomContainer oneMole) throws Exception{
		ECBasedBTransformer ecb = new ECBasedBTransformer(BioSystemName.HUMAN);
		ArrayList<Biotransformation> result = ecb.simulateECBasedPhaseIMetabolismChain(oneMole, true, true, 1, 0.0);
		return result;
	}
}
