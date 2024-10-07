/**
 * This class implements the class of metabolic enzymes.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.biomolecule;

import java.util.ArrayList;

import org.openscience.cdk.silent.SilentChemObjectBuilder;

import ambit2.smarts.SMIRKSManager;
import biotransformer.transformation.MetabolicReaction;

public class Enzyme {

	private String						name;
	private String						acceptedName;
	private ArrayList<String>				uniProtIds; //specific to a BioSystem
	private String 							primaryUniProtId;
	private ArrayList<MetabolicReaction>	reactionSet;
	private ArrayList<String> 				cellularLocations;
	private String							description;

	public Enzyme(String enzName) {
		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
		name = enzName;
//		reactionSet = generateReactionObjects(name, smrkMan);
	}


	public Enzyme(String enzName, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions) {	
//		Enzyme(enzName, null ,uniprot_ids, cellularLocations, reactions);
		this.name = enzName;
		this.uniProtIds = uniprot_ids;
		this.primaryUniProtId = uniprot_ids.get(0);
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		
	}

	
	
	public Enzyme(String enzName, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions, String acceptedcName) {	
//		Enzyme(enzName, null ,uniprot_ids, cellularLocations, reactions);
		this.name = enzName;
		this.acceptedName = acceptedcName;
		this.uniProtIds = uniprot_ids;
		this.primaryUniProtId = uniprot_ids.get(0);
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;

	}	
	public Enzyme(String enzName, String description, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions) {
		this.name = enzName;
		
		if(!(this.uniProtIds == null || this.uniProtIds.isEmpty())){
			this.uniProtIds = uniprot_ids;
			this.primaryUniProtId = uniprot_ids.get(0);
		}
		
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		this.description = description;
	}
	
	public Enzyme(String enzName, String description, ArrayList<String> uniprot_ids, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions, String systematicName) {
		this.name = enzName;
		this.acceptedName = acceptedName;
		
		if(!(this.uniProtIds == null || this.uniProtIds.isEmpty())){
			this.uniProtIds = uniprot_ids;
			this.primaryUniProtId = uniprot_ids.get(0);
		}
		
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		this.description = description;
	}
	
	public Enzyme(String enzName, String description, String uniprot_id, ArrayList<String> cellularLocations,
			ArrayList<MetabolicReaction> reactions) {
		this.name = enzName;
		this.uniProtIds = new ArrayList<String>();
		this.uniProtIds.add(uniprot_id);
		this.primaryUniProtId = uniprot_id;
		this.cellularLocations = cellularLocations;
		this.reactionSet = reactions;
		this.description = description;
		
	}

	public Enzyme(String enzName, ArrayList<MetabolicReaction>	reactionSet, ArrayList<String>cellLocations) {
		name = enzName;
		this.reactionSet = reactionSet;
		this.cellularLocations = cellLocations;

	}
	
	public Enzyme(String enzName, ArrayList<MetabolicReaction>	reactionSet, SMIRKSManager smrkMan) {
		name = enzName;
		this.reactionSet = reactionSet;
	}	
	
	public Enzyme(String enzName, ArrayList<MetabolicReaction>	reactionSet) {
		name = enzName;
		this.reactionSet = reactionSet;
	}
	


	
	public void addReaction(MetabolicReaction reaction) {
		reactionSet.add(reaction);
	}

	public String getName() {
		return name.toString();
	}

	public String getPrimaryUniprotId() {
		return primaryUniProtId;
	}

	public ArrayList<String> getUniprotIds() {
		return uniProtIds;
	}
	
	public ArrayList<MetabolicReaction> getReactionSet() {
		return reactionSet;
	}

	public void setprimaryUniprotId(String uniprotID) {
		primaryUniProtId = uniprotID;
	}
	
	public void setUniprotIds(ArrayList<String> uniprotIds) {
		uniProtIds =  uniprotIds;
	}


	public ArrayList<String> getReactionsNames() {
		ArrayList<String> names = new ArrayList<String>();
		for (int i = 0; i < reactionSet.size(); i++) {
			names.add(reactionSet.get(i).getReactionName());
		}
		return names;
	}
	
	@Override
	public String toString(){
		String representation = "";
		representation += String.format("%-20s\t%-25s\n","Name:",this.name);
		representation += String.format("%-30s\t%-35s\n","Systematic name:",this.acceptedName);
		representation += String.format("%-20s\n","Reactions:");
		
		for(MetabolicReaction m : this.getReactionSet()) {
			representation += String.format("%-60s\n",m.name);
		}
		
		return representation;
	}
	/**
	 * This function is used in the rails app, so we can get rid of that massive hardcoded things
	 * @param type
	 * @return
	 */
//	public String[] getEnzymeNameList(String type) throws Exception{
//		if(type.equalsIgnoreCase("CYP450")){
//			
//		}
//		else if(type.equalsIgnoreCase("ECBased")){
//			
//		}
//		else if(type.equalsIgnoreCase("PhaseII")){
//			Phase2BTransformer p2 = new Phase2BTransformer();
//		}
//		else if(type.equalsIgnoreCase("HGutBased")){
//			
//		}
//		else throw new Exception("Cannot find enzymes for the input biotransformer type");
//	}

}
