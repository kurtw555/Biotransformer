package biotransformer.dbrelevant;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.Utilities;

public class RetrieveFromDBHelpingFunctions {
	private GetBiotransformationType gbtype = new GetBiotransformationType();
	public final HashMap<String, String> queryTable;
	private ArrayList<String> inChIKeyList_ignore;
	public String data_source;
	private boolean first14;// = true;
	private boolean checkOverflow;// = true;
	private boolean substituted;// = false;
	private int overflow_threshold;
	private String query_biotransformation_type;
	private SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	public RetrieveFromDBHelpingFunctions(String source, String query_biotransformation_type, boolean first14, boolean checkOverflow, boolean substituted, int overflow_threshold) throws Exception{
		this.data_source = source;
		this.first14 = first14;
		this.checkOverflow = checkOverflow;
		this.substituted = substituted;
		this.overflow_threshold = overflow_threshold;
		this.inChIKeyList_ignore = inChIKeyToIgnore(first14);
		this.query_biotransformation_type = query_biotransformation_type;
		if(source.equalsIgnoreCase("pathbank")){
			this.queryTable = queryMap_Pathbank();
		}
		else {
			this.queryTable =  queryMap_hmdb();
		}
	}
	/**
	 * This function will retrieve metabolites for the query compound from Pathbank
	 * Here we add some extra restrctions:
	 * If there are more than 50 metabolites retrieved for the precursor, then we will return the first overflow_threshold metabolites with a comment for the reaction attribute saying that:
	 * there are too many metabolites, please run the jar file on your local machine to get the complete results.
	 * Please note that if checkOverflow is set as false, then overflow_threshold won't be used at all
	 * @param stmt
	 * @param inChIKey
	 * @param checkExported
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getBiotransformationInfo(Statement stmt, String inChIKey, boolean checkExported, boolean useSubstituition) throws Exception{		
		HashMap<String, HashMap<String, Object>> precursorEntry_map = getPrecusorEntryFromPathBank_byInChiKey(stmt, inChIKey, false, checkExported);
		String result_InChIKey = null;
		Integer precursorID;
		HashMap<Integer,ArrayList<Integer>> precursor_Reactions_map = new HashMap<>();
		if(!precursorEntry_map.isEmpty()){
			precursorID = (Integer) precursorEntry_map.get(inChIKey).get("Metabolite_ID");
			precursor_Reactions_map = getMetabolitesForReaction_usingMetaboliteID(stmt, precursorID);
		}		
		
		//If there are no metabolites for this compound, then use the first 14 charaters of the InChIKey to retrieve compounds with different(or no) stereo information
		if((precursorEntry_map.isEmpty() || precursor_Reactions_map.isEmpty()) && useSubstituition){
			System.out.println("Cannot find a record using the query InChIKey, hence we will try to find another matched record using the first 14 characters");
			//this.substituted = true;
			String inChIKey_14 = inChIKey.split("-")[0];
			precursorEntry_map = getPrecusorEntryFromPathBank_byInChiKey(stmt, inChIKey_14, true, checkExported);
			for(String inChIKey_alt : precursorEntry_map.keySet()){
				HashMap<String, Object> temp_map = precursorEntry_map.get(inChIKey_alt);
				Integer precursorID_temp = (Integer) temp_map.get("Metabolite_ID");
				HashMap<Integer,ArrayList<Integer>> precursor_Reactions_map_temp = getMetabolitesForReaction_usingMetaboliteID(stmt, precursorID_temp);	
				if(!precursor_Reactions_map_temp.isEmpty()){
					if(precursor_Reactions_map_temp.size() >  precursor_Reactions_map.size()){
						result_InChIKey = inChIKey_alt;
						precursor_Reactions_map = precursor_Reactions_map_temp;
						precursorID = precursorID_temp;
					}
				}
			}			
		}
		else result_InChIKey = inChIKey;
		if(result_InChIKey == null) return new ArrayList<Biotransformation>();
		HashMap<String, Object> precursorEntry_propertyMap =  precursorEntry_map.get(result_InChIKey);
		ArrayList<Biotransformation>  resultBiotransformationList = new ArrayList<>();
		for(Integer reactionID : precursor_Reactions_map.keySet()){
			if(this.checkOverflow){
				if(resultBiotransformationList.size() > this.overflow_threshold){
					markWithOverflowComment(resultBiotransformationList);
					break;
				}
			}
			ArrayList<Integer> metaboliteID_list = precursor_Reactions_map.get(reactionID);
			String[] enzyme_reaction = getEnzyme_usingReactionID(stmt, reactionID);
			String enzyme = enzyme_reaction[0];
			String reactionName = enzyme_reaction[1];
			if(enzyme.equals("N/A") || metaboliteID_list.isEmpty()) continue;
			ArrayList<Biotransformation> bt_temp_list = new ArrayList<>();
			boolean discardReaction = false;
			for(int i = 0; i < metaboliteID_list.size(); i++){
				int metaboliteID = metaboliteID_list.get(i);
				try{
					HashMap<String, Object> metabolitesEntry_propertyMap = getMetaboliteEntryFromHMDB_byMetaboliteID(stmt, metaboliteID, checkExported, reactionID);
					if(metabolitesEntry_propertyMap.isEmpty()) continue;
					//System.out.println((String)metabolitesEntry_propertyMap.get("moldb_smiles"));
					Biotransformation bt = generationBiotransformation(precursorEntry_propertyMap, metabolitesEntry_propertyMap, enzyme, reactionName);					
					//resultBiotransformationList.add(bt);
					//Should not discard the whole reactions if a metabolite is invalid. We should discard the invalid metabolite instead.
					if(bt == null){
						continue;
					}
					else{
						bt_temp_list.add(bt);
					}
				}catch(Exception e){
					//e.printStackTrace();
					System.out.println("Cannot process metabolite with CompoundID: " + metaboliteID);
					continue;
				}
			}
			if(!discardReaction) resultBiotransformationList.addAll(bt_temp_list);
			//System.out.println(String.format("ReactionID: %d, Enzyme: %s", reactionID, enzyme));
			//Get reaction information by reactionID. Then get Metabolite property Information by metaboliteID
			//Biotransformation bt = new Biotransformation(substrates, reactionType, enzymeName, products);
		}
		//if(this.substituted){
		if(useSubstituition) {
			markWithSubsitutionComment(resultBiotransformationList);
		}
		resultBiotransformationList = Utilities.selectUniqueBiotransformations(resultBiotransformationList);

		//String outputFileName = "E:/BiotransformerDependencies/Biotransformer3.0/temp_moloutput/Phosphoenolpyruvate_result.sdf";
		//ecb.saveBioTransformationProductsToSdf(resultBiotransformationList, outputFileName, new LinkedHashMap<String, MetabolicReaction>(), false);
		return resultBiotransformationList;
			
	}
	/**
	 * This function will form the biotransformation using the precursoor and the metabolite information
	 * If it returns null, then the whole reaction should be ignored.
	 * @param precursorMap
	 * @param metaboliteMap
	 * @param enzyme
	 * @return
	 * @throws Exception
	 */
	public Biotransformation generationBiotransformation(HashMap<String, Object> precursorMap, HashMap<String, Object> metaboliteMap, String enzyme, String reactionName) throws Exception{
		IAtomContainerSet precursorSet = generateAtomContainerObject(precursorMap);
		IAtomContainerSet productSet = generateAtomContainerObject(metaboliteMap);
		int precusorAtomCount = 0;
		int productAtomCount = 0;
		//If the metabolite has 21 or more atoms than the precursor,  then we think the precursor plays a very minor role in the reaction, and this reaction should be ignored.
		for(int i = 0; i < precursorSet.getAtomContainerCount(); i++){
			for(int j = 0; j < precursorSet.getAtomContainer(i).getAtomCount(); j++){
				if(!precursorSet.getAtomContainer(i).getAtom(j).getSymbol().equalsIgnoreCase("H")) precusorAtomCount++;
			}
		}
		for(int i = 0; i < productSet.getAtomContainerCount(); i++){
			if(!Utilities.isValidMetabolte(productSet.getAtomContainer(i))) return null;
			for(int j = 0; j < productSet.getAtomContainer(i).getAtomCount(); j++){
				if(!productSet.getAtomContainer(i).getAtom(j).getSymbol().equalsIgnoreCase("H")) productAtomCount++;
			}
			
		}
//		if(precusorAtomCount + 21 < productAtomCount){
//			return null;
//		}
		String reactionType = reactionName;
		ArrayList<String> enzymeList = new ArrayList<>();
		enzymeList.add(enzyme);
		Biotransformation biotransformation = new Biotransformation(precursorSet, reactionType, enzymeList, productSet, BioSystemName.HUMAN);
		return biotransformation;
 	}
	
	public IAtomContainerSet generateAtomContainerObject(HashMap<String, Object> moleculeMap) throws Exception{
		IAtomContainerSet  result = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		Integer metaboliteID = (Integer) moleculeMap.get("Metabolite_ID");
		String moldb_smiles = (String) moleculeMap.get("moldb_smiles");
		Double moldb_average_mass = (Double) moleculeMap.get("moldb_average_mass");
		Double moldb_mono_mass = (Double) moleculeMap.get("moldb_mono_mass");
		String moldb_formula = (String) moleculeMap.get("moldb_formula");
		Double moldb_alogps_logp = (Double) moleculeMap.get("moldb_alogps_logp");
		String moldb_inchikey = (String) moleculeMap.get("moldb_inchikey");
		String moldb_inchi = (String) moleculeMap.get("moldb_inchi");
		Integer reaction_id = (Integer) moleculeMap.get("Reaction ID");
		IAtomContainer structure = this.sp.parseSmiles(moldb_smiles);
		AtomContainerManipulator.suppressHydrogens(structure);
		
		structure.setProperty("InChI", moldb_inchi);
		structure.setProperty("InChIKey", moldb_inchikey);
		structure.setProperty("formula", moldb_formula);
		structure.setProperty("SMILES", moldb_smiles);
		structure.setProperty("Major Isotope Mass", String.valueOf(moldb_mono_mass));
		//structure.setProperty("RetrievedFromDB", true);
		if(reaction_id == null) structure.setProperty("Reaction ID", "N/A");
		else structure.setProperty("Reaction ID", String.valueOf(reaction_id));
		if(moldb_alogps_logp!=null){
			structure.setProperty("ALogP", moldb_alogps_logp.toString());
		}
		else{
			ALOGPDescriptor aLogpDescriptor = new ALOGPDescriptor(); 
			IDescriptorResult alogp = aLogpDescriptor.calculate(structure).getValue();
			structure.setProperty("ALogP", alogp.toString().split(",")[0]);
		}
		result.addAtomContainer(structure);
		return  result;
	}
	/**
	 * This function will get the enzyme name by reaciton ID
	 * @param stmt
	 * @param reactionID
	 * @return
	 * @throws Exception
	 */
	public String[] getEnzyme_usingReactionID(Statement stmt, Integer reactionID) throws Exception{
		String[] enzyme_reaction = new String[2];
		String query = String.format(this.queryTable.get("GetEnzyme_FromProteinTable_ByReactionID"), reactionID);
		ResultSet rs = stmt.executeQuery(query);
		if(rs.next()){
			String enzyme = rs.getString("name");
			String reactionName = rs.getString("reaction_name");
			if(enzyme == null || enzyme.equals("")) {
				enzyme = "N/A";
			}
			if(reactionName == null || reactionName.equals("")) {
				reactionName = "N/A";
			}
			enzyme_reaction[0] = enzyme;
			enzyme_reaction[1] = reactionName;
		}
		else{
			enzyme_reaction[0] = "N/A";
			enzyme_reaction[1] = "N/A";
			return enzyme_reaction;
			//throw new Exception("cannot find the enzyme for reacitonID: " + reactionID);
		}
		rs.close();
		return enzyme_reaction;
	}
	/**
	 * This function will retrieve all the metabolites stored in HMDB for the query metaboliteID_precursor
	 * It returns HashMap<ReactionID, MetaboliteID>
	 * It only considers those reactions that are not predicted by Biotransformer (predicted_in_hmdb) and discards those predicted ones.
	 * @param stmt
	 * @param metaboliteID_precursor
	 * @return
	 * @throws Exception
	 */
	public HashMap<Integer,ArrayList<Integer>> getMetabolitesForReaction_usingMetaboliteID(Statement stmt, Integer metaboliteID_precursor) throws Exception{
		String query = String.format(this.queryTable.get("GetReactionElement_FromReactonElementTable_ByMetaboliteID"), metaboliteID_precursor);
		ResultSet rs = stmt.executeQuery(query);
		ArrayList<Integer> validReactionID_list = new ArrayList<>();
		//Collect all reaction_ids for the reactions. We don't want to open multiple sql statments
		while(rs.next()){
			Integer reaction_id = rs.getInt("reaction_id");
			if(!validReactionID_list.contains(reaction_id)){				
				validReactionID_list.add(reaction_id);		
				//System.out.println("Reactoin_ID: " + reaction_id + "; Is predicted by biotransforemr: " + isNotPredicted);
				
			}							
		}
		rs.close();
		for(int i = 0; i < validReactionID_list.size(); i++) {
			if(!this.gbtype.isQueriedBiotransformation(stmt, validReactionID_list.get(i), this.query_biotransformation_type)) {
				validReactionID_list.remove(i);
				i = i-1;
			}
		}
		if(this.data_source.equals("hmdb")){
			for(int i = 0; i < validReactionID_list.size(); i++){	
				//If the reaction is predicted by BioTransformer, then remove it from the list
				if(!isNotPredicted_From_Reaction_Table_ByReactionID(stmt,validReactionID_list.get(i))){
					validReactionID_list.remove(i);
					i = i-1;
				}
			}
		}
		
		HashMap<Integer,ArrayList<Integer>> reactionMetabolite_map = new HashMap<>();
		for(int i = 0; i < validReactionID_list.size(); i++){
			boolean isValid = true;
			//System.out.println("Reaction ID: " + validReactionID_list.get(i));
			ArrayList<Integer> metaboliteID_oneReaction = new ArrayList<>();
			query = String.format(this.queryTable.get("GetReactionElement_FromReactonElementTable_ByReactionID"),validReactionID_list.get(i));
			//Add all metabolites in one reaction
			rs = stmt.executeQuery(query);
			while(rs.next()){
				Integer metabolite_id = rs.getInt("element_id");
				if(metabolite_id == null || metabolite_id == 0) continue;
				if(this.data_source.equalsIgnoreCase("pathbank")){
					String element_type = rs.getString("element_type");
					if(!element_type.equals("Compound")){
						isValid = false;
						break;
					}
				}
				metaboliteID_oneReaction.add(metabolite_id);
			}
			if(!isIncludedSet(reactionMetabolite_map.values(), metaboliteID_oneReaction) && isValid){
				reactionMetabolite_map.put(validReactionID_list.get(i), metaboliteID_oneReaction);
			}
			
			rs.close();
		}		
		return reactionMetabolite_map;
	}	
	
	/**
	 * This function will get the molecule(metabolite) information for a compound by its metaboliteID
	 * This function will only be used when retrieving the information of a metabolite in one reaction which is stored in HMDBs
	 * @param stmt
	 * @param metaboliteID
	 * @return
	 * @throws Exception
	 */
	public HashMap<String,Object> getMetaboliteEntryFromHMDB_byMetaboliteID(Statement stmt, Integer metaboliteID, boolean checkExported, Integer reactionID) throws Exception{
		//System.out.println(metaboliteID);
		String query = String.format(this.queryTable.get("GetMetabolite_FromMetaboliteTatble_ByMetaboliteID"), metaboliteID);
		ResultSet rs = stmt.executeQuery(query);
		rs.next();
		HashMap<String, Object> resultMap_oneInChI = new HashMap<>();
		if(this.data_source.equalsIgnoreCase("hmdb")){
			boolean export_to_hmdb = rs.getBoolean("export_to_hmdb");
			if(!export_to_hmdb && checkExported) return resultMap_oneInChI;
		}
		String moldb_smiles = rs.getString("moldb_smiles");
		Double moldb_average_mass = rs.getDouble("moldb_average_mass");
		Double moldb_mono_mass = rs.getDouble("moldb_mono_mass");
		String moldb_formula = rs.getString("moldb_formula");
		Double moldb_alogps_logp = rs.getDouble("moldb_alogps_logs");
		String moldb_inchikey = rs.getString("moldb_inchikey");
		String moldb_inchi = rs.getString("moldb_inchi");
		resultMap_oneInChI.put("Metabolite_ID", metaboliteID);
		resultMap_oneInChI.put("moldb_smiles", moldb_smiles);
		resultMap_oneInChI.put("moldb_average_mass", moldb_average_mass);
		resultMap_oneInChI.put("moldb_mono_mass", moldb_mono_mass);
		resultMap_oneInChI.put("moldb_formula", moldb_formula);
		resultMap_oneInChI.put("moldb_alogps_logp", moldb_alogps_logp);
		resultMap_oneInChI.put("moldb_inchikey", moldb_inchikey);	
		resultMap_oneInChI.put("moldb_inchi", moldb_inchi);	
		resultMap_oneInChI.put("Reaction ID", reactionID);
		if(this.inChIKeyList_ignore.contains(moldb_inchikey.split("-")[0])) return new HashMap<String, Object>();
		IAtomContainer molecule = this.sp.parseSmiles(moldb_smiles);
		int countNonHydrogen = 0;
		for(int i = 0; i < molecule.getAtomCount(); i++){
			if(!molecule.getAtom(i).getSymbol().equalsIgnoreCase("H")) countNonHydrogen++;
			if(countNonHydrogen > 1) break;
		}
		if(countNonHydrogen<=1) return new HashMap<String, Object>();
		//System.out.println(moldb_inchi);
		return resultMap_oneInChI;
	}
	/**
	 * This function will get the molecule information using the inChiKey
	 * Return HashMap<InChiKey, HashMap<Property, value>>
	 * @param stmt
	 * @param inChIKey
	 * @return
	 * @throws Exception
	 */
	public HashMap<String, HashMap<String,Object>> getPrecusorEntryFromPathBank_byInChiKey(Statement stmt, String inChIKey, boolean usefirst14, boolean checkExported) throws Exception{
		HashMap<String, HashMap<String, Object>> resultMap = new HashMap<>();				
		String query;		
		if(!usefirst14){
			query = String.format(this.queryTable.get("GetSubstrate_FromMetaboliteTable_ByInChIKey"), inChIKey);			
		}
		else{
			query = String.format(this.queryTable.get("GetSubstrate_FromMetaboliteTable_ByInChIKey_14"), inChIKey);
		}
		ResultSet rs = stmt.executeQuery(query);
		while(rs.next()){
			HashMap<String, Object> resultMap_oneInChI = new HashMap<>();
			Integer metaboliteID = rs.getInt("id");
			if(this.data_source.equalsIgnoreCase("hmdb")){
				boolean export_to_hmdb = rs.getBoolean("export_to_hmdb");
				if(!export_to_hmdb && checkExported) continue;
			}
			String moldb_smiles = rs.getString("moldb_smiles");
			Double moldb_average_mass = rs.getDouble("moldb_average_mass");
			Double moldb_mono_mass = rs.getDouble("moldb_mono_mass");
			String moldb_formula = rs.getString("moldb_formula");
			Double moldb_alogps_logp = rs.getDouble("moldb_alogps_logs");
			String moldb_inchikey = rs.getString("moldb_inchikey");
			String moldb_inchi = rs.getString("moldb_inchi");
			//Close this result set every time after it is used			
			resultMap_oneInChI.put("Metabolite_ID", metaboliteID);
			resultMap_oneInChI.put("moldb_smiles", moldb_smiles);
			resultMap_oneInChI.put("moldb_average_mass", moldb_average_mass);
			resultMap_oneInChI.put("moldb_mono_mass", moldb_mono_mass);
			resultMap_oneInChI.put("moldb_formula", moldb_formula);
			resultMap_oneInChI.put("moldb_alogps_logp", moldb_alogps_logp);
			resultMap_oneInChI.put("moldb_inchikey", moldb_inchikey);
			resultMap_oneInChI.put("moldb_inchi", moldb_inchi);
			//System.out.println("Processing MetaboliteID: " + metaboliteID);
			//If two entries have the same InChIKey but different metaboliteID, then take the smaller one
			if(resultMap.containsKey(moldb_inchikey)){
				Integer tempID = (Integer) resultMap.get(moldb_inchikey).get("Metabolite_ID");
				if(metaboliteID < tempID) resultMap.put(moldb_inchikey, resultMap_oneInChI);
			}
			else{
				resultMap.put(moldb_inchikey, resultMap_oneInChI);
			}
		}
		rs.close();
		return resultMap;
		
	}
	/**
	 * This function will use a query to check if the retrieved reaction is predicted by BioTransformer
	 * Deprecated in RetrieveFromPathBank module
	 * @param stmt
	 * @return
	 * @throws Exception
	 */
	public boolean isNotPredicted_From_Reaction_Table_ByReactionID(Statement stmt, Integer reactionID) throws Exception{
		String query = String.format(this.queryTable.get("GetReaction_FromReactionTable_ByReactionID"), reactionID);
		ResultSet rs = stmt.executeQuery(query);
		if(rs.getFetchSize() > 1){
			throw new Exception("There are more than one records matching the query ID in the isNotPredicted_From_Reaction_Table_ByReactionID function");
		}
		rs.next();
//		System.out.println("Reaction_ID: " + rs.getInt("id"));
//        System.out.println("altext: " + rs.getString("altext"));
		boolean result = rs.getBoolean("predicted_in_hmdb");		
		//Close this result set every time after it is used
		rs.close();
		return !result;
	}
	/**
	 * This function checks if the set_one has been included in the metaboliteSets.
	 * @param metaboliteSets
	 * @param set_one
	 * @return
	 */
	public boolean isIncludedSet(Collection<ArrayList<Integer>> metaboliteSets, ArrayList<Integer> set_one){
		Iterator<ArrayList<Integer>> iterate_metaboliteSets = metaboliteSets.iterator();
		while(iterate_metaboliteSets.hasNext()){
			ArrayList<Integer> temp = iterate_metaboliteSets.next();
			if(areSameMetaboliteSet(temp, set_one)) return true;
		}
		return false;
	}
	/**
	 * This function checks if two sets of metabolites are identical by their metaboliteIDs
	 * @param set_one
	 * @param set_two
	 * @return
	 */
	public boolean areSameMetaboliteSet(ArrayList<Integer> set_one, ArrayList<Integer> set_two){
		if(set_one.size() != set_two.size()) return false;
		for(int i = 0; i < set_one.size(); i++){
			if(!set_two.contains(set_one.get(i))) return false;
		}
		return true;
	}
	/**
	 * This function creates a HashMap that contains all the <QueryName, QueryText> maps.
	 * The QueryName and QueryText are predefined are used to retrieve data from the PAthbank database.
	 * @return 
	 */
	public HashMap<String, String> queryMap_Pathbank(){
		HashMap<String, String> resultMap = new HashMap<>();		
		resultMap.put("GetReaction_FromReactionTable_ByReactionID", "select * from pathbank.reactions where id in (select reaction_id from pathbank.reaction_elements where element_type = 'Compound' and type = 'ReactionLeftElement' and element_id = '%d')");
		resultMap.put("GetSubstrate_FromMetaboliteTable_ByInChIKey_14", "select id, moldb_inchikey,moldb_alogps_logs, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula FROM pathbank.compounds where moldb_inchikey like '%s%%'");
		resultMap.put("GetSubstrate_FromMetaboliteTable_ByInChIKey", "select id, moldb_inchikey,moldb_alogps_logs, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula FROM pathbank.compounds where moldb_inchikey = '%s'");
		resultMap.put("GetMetabolite_FromMetaboliteTatble_ByMetaboliteID", "select id, moldb_inchikey,moldb_alogps_logs, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula from pathbank.compounds where id = %d");
		resultMap.put("GetReactionElement_FromReactonElementTable_ByMetaboliteID", "select reaction_id from pathbank.reaction_elements where element_type = 'Compound' and type = 'ReactionLeftElement' and element_id = %d");
		resultMap.put("GetReactionElement_FromReactonElementTable_ByReactionID", "select element_id, element_type from pathbank.reaction_elements where reaction_id = %d and type = 'ReactionRightElement'");	
		resultMap.put("GetEnzyme_FromProteinTable_ByReactionID", "select * from pathbank.protein_complexes where id in (select protein_complex_id from pathbank.reaction_enzymes where reaction_id = %d)");
		return resultMap;
	}
	/**
	 * This function creates a HashMap that contains all the <QueryName, QueryText> maps.
	 * The QueryName and QueryText are predefined are used to retrieve data from the HMDB database.
	 * @return 
	 */
	public HashMap<String, String> queryMap_hmdb(){
		HashMap<String, String> resultMap = new HashMap<>();
		resultMap.put("GetReaction_FromReactionTable_ByReactionID", "select * from hmdb.reactions where id = '%d'");
		resultMap.put("GetSubstrate_FromMetaboliteTable_ByInChIKey_14", "select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp as moldb_alogps_logs from hmdb.metabolites where moldb_inchikey like '%s%%'");
		resultMap.put("GetSubstrate_FromMetaboliteTable_ByInChIKey", "select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp as moldb_alogps_logs from hmdb.metabolites where moldb_inchikey = '%s'");
		resultMap.put("GetMetabolite_FromMetaboliteTatble_ByMetaboliteID", "select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp as moldb_alogps_logs from hmdb.metabolites where id = %d");
		resultMap.put("GetReactionElement_FromReactonElementTable_ByMetaboliteID", "select reaction_id from hmdb.reaction_elements where metabolite_id = '%d' AND type = 'LeftReactionElement'");
		resultMap.put("GetReactionElement_FromReactonElementTable_ByReactionID", "select metabolite_id as element_id from hmdb.reaction_elements where reaction_id = %d AND type = 'RightReactionElement'");	
		resultMap.put("GetEnzyme_FromProteinTable_ByReactionID", "select name, reaction_name from reaction_names join proteins on proteins.enzyme_classification = reaction_names.enzyme_classification where proteins.id in (select protein_id from reaction_enzymes where reaction_id = %d)");
		resultMap.put("GetMetabolitesContainsElements", "select * from hmdb.metabolites, hmdb.ontologies where export_to_hmdb = 1 AND moldb_smiles like '%s%%' AND ontologies.metabolite_id = hmdb.metabolites.id AND ontologies.origin = 'Endogenous' AND detected = 1");
		return resultMap;
	}

	/**
	 * This creates the arrayList that contains all the InChIKeys of those CoEnzymeA-ish compounds that should not be considered.
	 * If boolean first14Characters is true, then it will return the ArrayList<String> that contains the first 14 characters of the InChIKeys
	 * @return
	 */
	public ArrayList<String> inChIKeyToIgnore(boolean first14Characters){
		ArrayList<String> resultList = new ArrayList<>();
		resultList.add("ACFIXJIJDZMPPO-NNYOXOHSSA-N");
		resultList.add("ACTIUHUUMQJHFO-UPTCCGCDSA-N");
		resultList.add("AZQWKYJCGOJGHM-UHFFFAOYSA-N");
		resultList.add("BAWFJGJZGIEFAR-NNYOXOHSSA-N");
		resultList.add("BJSUUWFQAMLNKU-OKXKTURISA-N");
		resultList.add("BOPGDPNILDQYTO-NNYOXOHSSA-N");
		resultList.add("CIKGWCTVFSRMJU-KVQBGUIXSA-N");
		resultList.add("IVOMOUWHDPKRLL-KQYNXXCUSA-N");
		resultList.add("NBIIXXVUZAFLBC-UHFFFAOYSA-N");
		resultList.add("QGWNDRXFNXRZMB-UUOKFMHZSA-N");
		resultList.add("RGJOEKWQDUBAIZ-IBOSZNHHSA-N");
		resultList.add("RNUCUWWMTTWKAH-JLHYYAGUSA-N");
		resultList.add("SENPVEZBRZQVST-HISDBWNOSA-O");
		resultList.add("UDMBCSSLTHHNCD-KQYNXXCUSA-N");
		resultList.add("UNXRWKVEANCORM-UHFFFAOYSA-N");
		resultList.add("VWWQXMAJTJZDQX-UYBVJOGSSA-N");
		resultList.add("WHTCPDAXWFLDIH-KQYNXXCUSA-N");
		resultList.add("XJLXINKUBYWONI-NNYOXOHSSA-N");
		resultList.add("XKMLYUALXHKNFT-UUOKFMHZSA-N");
		resultList.add("XPPKVPWEQAFLFU-UHFFFAOYSA-N");
		resultList.add("XTWYTFMLZFPYCI-KQYNXXCUSA-N");
		resultList.add("YPZRHBJKEMOYQH-UYBVJOGSSA-N");
		resultList.add("YPZRHBJKEMOYQH-UYBVJOGSSA-N");
		resultList.add("ZKHQWZAMYRWXGA-KQYNXXCUSA-N");
		ArrayList<String> first14List = new ArrayList<>();
		for(int i = 0; i < resultList.size(); i++){
			if(!first14List.contains(resultList.get(i).split("-")[0])) first14List.add(resultList.get(i).split("-")[0]);
		}
		if(first14Characters) return first14List;
		return resultList;
	}
	
	/**
	 * This function takes the ArrayList<Biotransformation> that is 30 or longer for a precursor 
	 * @param resultList
	 * @return
	 * @throws Exception
	 */
	ArrayList<Biotransformation> markWithOverflowComment(ArrayList<Biotransformation> resultList) throws Exception{
		for(int i = 0; i < resultList.size(); i++){
			Biotransformation bt = resultList.get(i);
			bt.setReactionType(bt.getReactionType() + ";" + "*Overflow");						
		}
		return resultList;
	}
	/**
	 * This function assigns a comment that indicates the precusor is substituted by another record in HMDB for that Biotransformation 
	 * @param resultList
	 * @return
	 * @throws Exception
	 */
	ArrayList<Biotransformation> markWithSubsitutionComment(ArrayList<Biotransformation> resultList) throws Exception{
		for(int i = 0; i < resultList.size(); i++){
			Biotransformation bt = resultList.get(i);
			bt.setReactionType(bt.getReactionType() + ";" + "*SubstitutedInChIKey");
		}
		return resultList;
	}
	
	public ArrayList<String> getInChIKeyList_ignore(){
		return this.inChIKeyList_ignore;
	}
	
	
}
