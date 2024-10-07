package biotransformer.utils;

import java.util.ArrayList;

import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;

public class GetEnzymeList {
	boolean useDB = false;
	boolean useSubstitution = false;
	public static enum BioType {
		CYP450, ECBased, PhaseII, HGutBased, Environmental};
	public static void main(String[] args) throws Exception{
		GetEnzymeList gel = new GetEnzymeList();
		String[] result = gel.getEnzymeList(BioType.CYP450.toString());
		for(int i = 0; i < result.length; i++){
			System.out.println(result[i]);
		}
	}
	/**
	 * This function is used to replace the hardcoded section in the biotransformerrails code.
	 * It generates the enzymeLists for the query biotransformer type automatically 
	 * @param type
	 * @return
	 * @throws Exception
	 */
	public String[] getEnzymeList(String type) throws Exception{
		ArrayList<Enzyme> enzymeArray;
		if(type.equalsIgnoreCase("CYP450")){
			Cyp450BTransformer cyp450 = new Cyp450BTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
			enzymeArray = cyp450.enzymesList;
		}
		else if(type.equalsIgnoreCase("ECBased")){
			ECBasedBTransformer ecb = new ECBasedBTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
			enzymeArray = ecb.enzymesList;
		}
		else if(type.equalsIgnoreCase("PhaseII")){
			Phase2BTransformer p2 = new Phase2BTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
			enzymeArray = p2.enzymesList;
		}
		else if(type.equalsIgnoreCase("HGutBased")){
			HGutBTransformer hgut = new HGutBTransformer(useDB, useSubstitution);
			enzymeArray = hgut.enzymesList;
		}
		else if(type.equalsIgnoreCase("Environmental")){
			EnvMicroBTransformer env = new EnvMicroBTransformer();
			enzymeArray = env.enzymesList;
		}
		else throw new Exception("Cannot find enzymes for the input biotransformer type");
		
		String[] enzymeList = new String[enzymeArray.size()];
		for(int i = 0; i < enzymeArray.size(); i++){
			enzymeList[i] = enzymeArray.get(i).getName();
		}
		return enzymeList;
	}
	/**
	 * This is an alternative function that takes the BioType as input and converts the enum to String internally
	 * @param typeEnum
	 * @return
	 * @throws Exception
	 */
	public String[] getEnzymeList(BioType typeEnum) throws Exception{
		ArrayList<Enzyme> enzymeArray;
		String type = typeEnum.toString();
		if(type.equalsIgnoreCase("CYP450")){
			Cyp450BTransformer cyp450 = new Cyp450BTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
			enzymeArray = cyp450.enzymesList;
		}
		else if(type.equalsIgnoreCase("ECBased")){
			ECBasedBTransformer ecb = new ECBasedBTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
			enzymeArray = ecb.enzymesList;
		}
		else if(type.equalsIgnoreCase("PhaseII")){
			Phase2BTransformer p2 = new Phase2BTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
			enzymeArray = p2.enzymesList;
		}
		else if(type.equalsIgnoreCase("HGutBased")){
			HGutBTransformer hgut = new HGutBTransformer(useDB, useSubstitution);
			enzymeArray = hgut.enzymesList;
		}
		else if(type.equalsIgnoreCase("Environmental")){
			EnvMicroBTransformer env = new EnvMicroBTransformer();
			enzymeArray = env.enzymesList;
		}
		else throw new Exception("Cannot find enzymes for the input biotransformer type");
		
		String[] enzymeList = new String[enzymeArray.size()];
		for(int i = 0; i < enzymeArray.size(); i++){
			enzymeList[i] = enzymeArray.get(i).getName();
		}
		return enzymeList;
	}
}
