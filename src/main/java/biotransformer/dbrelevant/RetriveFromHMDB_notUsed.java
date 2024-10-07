package biotransformer.dbrelevant;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.Utilities;

/**
 * Because most reactions/enzymes in HMDB are EC based and the PhaseI/PhaseII/GutMicrobial information are not given, this RetriveFromHMDB class is used as a module in ECBased biotransformation prediction
 * @author Tian
 *
 */
public class RetriveFromHMDB_notUsed {
	private final HashMap<String, String> queryTable;
	private final String DB_URL = "jdbc:mysql://162.243.63.130:3306/hmdb";
	private final String USER = "readonly";
	private final String PASS = "tardis";
	private SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	
	public static void main(String[] args) throws Exception{
		RetriveFromHMDB_notUsed rfdb = new RetriveFromHMDB_notUsed();
		boolean checkExported = true;
		//String InChIKey = "MUMGGOZAMZWBJJ-UHFFFAOYSA-N";
		//String InChIKey = "PFTAWBLQPZVEMU-UKRRQHHQSA-N";
		//String InChIKey = "WQZGKKKJIJFFOK-VFUOTHLCSA-N";
		String InChIKey = "RGJOEKWQDUBAIZ-IBOSZNHHSA-N";
		ArrayList<Biotransformation> result = rfdb.getPredictedBiotransformation(InChIKey, checkExported);
		System.out.println(result.size());
	}
	
	public RetriveFromHMDB_notUsed(){
		this.queryTable = queryMap();
	}
	/**
	 * This is the function that should be called in the ECBased metabolism prediction
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getBiotransformationRetrievedFromHMDB(IAtomContainer molecule, boolean checkExported) throws Exception{
		InChIGeneratorFactory inChIGenerator = InChIGeneratorFactory.getInstance();
		InChIGenerator inchiGen = inChIGenerator.getInChIGenerator(molecule);
		String inChiKey = inchiGen.getInchiKey();
		ArrayList<Biotransformation> results = getPredictedBiotransformation(inChiKey, checkExported);
		return results;
	}
	/**
	 * This function will be used in the ruby rails app and return ArrayList<Biotransformation>
	 * @param inChIKey
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getPredictedBiotransformation(String inChiKey, boolean checkExported) throws Exception{
		ArrayList<Biotransformation> resultArrayList = new ArrayList<>();
		try(
				Connection conn = DriverManager.getConnection(this.DB_URL, this.USER, this.PASS);
				Statement stmt = conn.createStatement();
			    ){		      

					ArrayList<Biotransformation> oneResult = getBiotransformationInfo(stmt, inChiKey, checkExported);
					resultArrayList.addAll(oneResult);
			      } catch (Exception e) {
					e.printStackTrace();
				}
		//System.out.println(resultArrayList.size());
		return resultArrayList;
	}
	public ArrayList<Biotransformation> getBiotransformationInfo(Statement stmt, String inChIKey, boolean checkExported) throws Exception{
		
		HashMap<String, HashMap<String, Object>> precursorEntry_map = getPrecusorEntryFromHMDB_byInChiKey(stmt, inChIKey, false, checkExported);
		String result_InChIKey = null;
		Integer precursorID = (Integer) precursorEntry_map.get(inChIKey).get("Metabolite_ID");
		HashMap<Integer,ArrayList<Integer>> precursor_Reactions_map = getReaction_usingMetaboliteID(stmt, precursorID);
		//If there are no metabolites for this compound, then use the first 14 charaters of the InChIKey to retrieve compounds with different(or no) stereo information
		if(precursor_Reactions_map.isEmpty()){
			String inChIKey_14 = inChIKey.split("-")[0];
			precursorEntry_map = getPrecusorEntryFromHMDB_byInChiKey(stmt, inChIKey_14, true, checkExported);
			for(String inChIKey_alt : precursorEntry_map.keySet()){
				HashMap<String, Object> temp_map = precursorEntry_map.get(inChIKey_alt);
				precursorID = (Integer) temp_map.get("Metabolite_ID");
				precursor_Reactions_map = getReaction_usingMetaboliteID(stmt, precursorID);	
				if(!precursor_Reactions_map.isEmpty()){
					result_InChIKey = inChIKey_alt;
					break;
				}
			}			
		}
		else result_InChIKey = inChIKey;
		if(result_InChIKey == null) return new ArrayList<Biotransformation>();
		HashMap<String, Object> precursorEntry_propertyMap =  precursorEntry_map.get(result_InChIKey);
		System.out.println("Final InChiKey used: " + precursorEntry_propertyMap.get("moldb_inchikey"));
		System.out.println("Final SMILES Strings used: " + precursorEntry_propertyMap.get("moldb_smiles"));
		ArrayList<Biotransformation>  resultBiotransformationList = new ArrayList<>();
		for(Integer reactionID : precursor_Reactions_map.keySet()){
			ArrayList<Integer> metaboliteID_list = precursor_Reactions_map.get(reactionID);
			String enzyme = getEnzyme_usingReactionID(stmt, reactionID);
			for(int i = 0; i < metaboliteID_list.size(); i++){
				int metaboliteID = metaboliteID_list.get(i);
				HashMap<String, Object> metabolitesEntry_propertyMap = getMetaboliteEntryFromHMDB_byMetaboliteID(stmt, metaboliteID, checkExported);
				//System.out.println((String)metabolitesEntry_propertyMap.get("moldb_smiles"));
				Biotransformation bt = generationBiotransformation(precursorEntry_propertyMap, metabolitesEntry_propertyMap, enzyme);
				resultBiotransformationList.add(bt);
			}
			//System.out.println(String.format("ReactionID: %d, Enzyme: %s", reactionID, enzyme));
			//Get reaction information by reactionID. Then get Metabolite property Information by metaboliteID
			//Biotransformation bt = new Biotransformation(substrates, reactionType, enzymeName, products);
		}
		ECBasedBTransformer ecb	= new ECBasedBTransformer(BioSystemName.HUMAN, true, true);
		String outputFileName = "E:/BiotransformerDependencies/Biotransformer3.0/temp_moloutput/testHMDBConnection_result.sdf";
		ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(resultBiotransformationList), outputFileName, new LinkedHashMap<String, MetabolicReaction>(), false);
		return resultBiotransformationList;
			
	}
	
	public Biotransformation generationBiotransformation(HashMap<String, Object> precursorMap, HashMap<String, Object> metaboliteMap, String enzyme) throws Exception{
		IAtomContainerSet precursorSet = generateAtomContainerObject(precursorMap);
		IAtomContainerSet productSet = generateAtomContainerObject(metaboliteMap);
		String temp = enzyme.substring(0, enzyme.indexOf("ase") + 3);
		String reactionType = temp.replace("tase", "tion").replace("nase", "nation").replace("dase", "dation").replace("rase", "ration");
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
		
		IAtomContainer structure = this.sp.parseSmiles(moldb_smiles);
		AtomContainerManipulator.suppressHydrogens(structure);
		
		structure.setProperty("InChI", moldb_inchi);
		structure.setProperty("InChIKey", moldb_inchikey);
		structure.setProperty("formula", moldb_formula);
		structure.setProperty("SMILES", moldb_smiles);
		structure.setProperty("Major Isotope Mass", String.valueOf(moldb_mono_mass));
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
	public String getEnzyme_usingReactionID(Statement stmt, Integer reactionID) throws Exception{
		String query = String.format(this.queryTable.get("GetEnzyme_FromProteinTable_ByReactionID"), reactionID);
		ResultSet rs = stmt.executeQuery(query);
		if(rs.next()){
			String enzyme = rs.getString("name");
			return enzyme;
		}
		else throw new Exception("cannot find the enzyme for reacitonID: " + reactionID);
	}
	/**
	 * This function will retrieve all the metabolites stored in HMDB for the query metaboliteID_precursor
	 * It only considers those reactions that are not predicted by Biotransformer (predicted_in_hmdb) and discards those predicted ones.
	 * @param stmt
	 * @param metaboliteID_precursor
	 * @return
	 * @throws Exception
	 */
	public HashMap<Integer,ArrayList<Integer>> getReaction_usingMetaboliteID(Statement stmt, Integer metaboliteID_precursor) throws Exception{
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
		//Remove the reactions that are predicted by BioTransformer
		for(int i = 0; i < validReactionID_list.size(); i++){	
			//If the reaction is predicted by BioTransformer, then remove it from the list
			if(!isNotPredicted_From_Reaction_Table_ByReactionID(stmt,validReactionID_list.get(i))){
				validReactionID_list.remove(i);
				i = i-1;
			}
		}
		HashMap<Integer,ArrayList<Integer>> reactionMetabolite_map = new HashMap<>();
		for(int i = 0; i < validReactionID_list.size(); i++){
			//System.out.println("Reaction ID: " + validReactionID_list.get(i));
			ArrayList<Integer> metaboliteID_oneReaction = new ArrayList<>();
			query = String.format(this.queryTable.get("GetReactionElement_FromReactonElementTable_ByReactionID"),validReactionID_list.get(i));
			//Add all metabolites in one reaction
			rs = stmt.executeQuery(query);
			while(rs.next()){
				Integer metabolite_id = rs.getInt("metabolite_id");
				if(metabolite_id == null || metabolite_id == 0) continue;
				metaboliteID_oneReaction.add(metabolite_id);
			}
			if(!isIncludedSet(reactionMetabolite_map.values(), metaboliteID_oneReaction)){
				reactionMetabolite_map.put(validReactionID_list.get(i), metaboliteID_oneReaction);
			}
			
			rs.close();
		}		
//		System.out.println("Total number of Reactions: " + validReactionID_list.size());
//		for(Integer reactionID : reactionMetabolite_map.keySet()){
//			System.out.println("Reaction ID: " + reactionID);
//			ArrayList<Integer> metaboliteList = reactionMetabolite_map.get(reactionID);
//			for(int i = 0; i < metaboliteList.size(); i++){
//				System.out.print("MetaboliteID: " + metaboliteList.get(i) + " ");
//			}
//			System.out.println();
//			System.out.println("-----------------------");
//		}
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
	public HashMap<String,Object> getMetaboliteEntryFromHMDB_byMetaboliteID(Statement stmt, Integer metaboliteID, boolean checkExported) throws Exception{
		//System.out.println(metaboliteID);
		String query = String.format(this.queryTable.get("GetMetabolite_FromMetaboliteTatble_ByMetaboliteID"), metaboliteID);
		ResultSet rs = stmt.executeQuery(query);
		rs.next();
		HashMap<String, Object> resultMap_oneInChI = new HashMap<>();
		boolean export_to_hmdb = rs.getBoolean("export_to_hmdb");
		if(!export_to_hmdb && checkExported) return resultMap_oneInChI;
		String moldb_smiles = rs.getString("moldb_smiles");
		Double moldb_average_mass = rs.getDouble("moldb_average_mass");
		Double moldb_mono_mass = rs.getDouble("moldb_mono_mass");
		String moldb_formula = rs.getString("moldb_formula");
		Double moldb_alogps_logp = rs.getDouble("moldb_alogps_logp");
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
	public HashMap<String, HashMap<String,Object>> getPrecusorEntryFromHMDB_byInChiKey(Statement stmt, String inChIKey, boolean usefirst14, boolean checkExported) throws Exception{
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
			boolean export_to_hmdb = rs.getBoolean("export_to_hmdb");
			if(!export_to_hmdb && checkExported) continue;
			String moldb_smiles = rs.getString("moldb_smiles");
			Double moldb_average_mass = rs.getDouble("moldb_average_mass");
			Double moldb_mono_mass = rs.getDouble("moldb_mono_mass");
			String moldb_formula = rs.getString("moldb_formula");
			Double moldb_alogps_logp = rs.getDouble("moldb_alogps_logp");
			String moldb_inchikey = rs.getString("moldb_inchikey");
			String moldb_inchi = rs.getString("moldb_inchi");
	//		System.out.println("Metabolite_ID: " + metaboliteID);
	//		System.out.println("moldb_smiles: " + moldb_smiles);		
	//		System.out.println("moldb_average_mass: " + moldb_average_mass);
	//		System.out.println("moldb_mono_mass: " + moldb_mono_mass);
	//		System.out.println("moldb_formula: " + moldb_formula);
	//		System.out.println("moldb_alogps_logp: " + moldb_alogps_logp);
	//		System.out.println("moldb_inchikey: " + moldb_inchikey);
	//		System.out.println();	
			//Close this result set every time after it is used			
			resultMap_oneInChI.put("Metabolite_ID", metaboliteID);
			resultMap_oneInChI.put("moldb_smiles", moldb_smiles);
			resultMap_oneInChI.put("moldb_average_mass", moldb_average_mass);
			resultMap_oneInChI.put("moldb_mono_mass", moldb_mono_mass);
			resultMap_oneInChI.put("moldb_formula", moldb_formula);
			resultMap_oneInChI.put("moldb_alogps_logp", moldb_alogps_logp);
			resultMap_oneInChI.put("moldb_inchikey", moldb_inchikey);
			resultMap_oneInChI.put("moldb_inchi", moldb_inchi);
			resultMap.put(moldb_inchikey, resultMap_oneInChI);
		}
		rs.close();
		return resultMap;
		
	}
	/**
	 * This function will use a query to check if the retrieved reaction is predicted by BioTransformer
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
	 * The QueryName and QueryText are predefined are used to retrieve data from the HMDB database.
	 * @return 
	 */
	public HashMap<String, String> queryMap(){
		HashMap<String, String> resultMap = new HashMap<>();
		resultMap.put("GetReaction_FromReactionTable_ByReactionID", "select * from hmdb.reactions where id = '%d'");
		resultMap.put("GetSubstrate_FromMetaboliteTable_ByInChIKey_14", "select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp from hmdb.metabolites where moldb_inchikey like '%s%%'");
		resultMap.put("GetSubstrate_FromMetaboliteTable_ByInChIKey", "select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp from hmdb.metabolites where moldb_inchikey = '%s'");
		resultMap.put("GetMetabolite_FromMetaboliteTatble_ByMetaboliteID", "select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp from hmdb.metabolites where id = %d");
		resultMap.put("GetReactionElement_FromReactonElementTable_ByMetaboliteID", "select reaction_id from hmdb.reaction_elements where metabolite_id = '%d' AND type = 'LeftReactionElement'");
		resultMap.put("GetReactionElement_FromReactonElementTable_ByReactionID", "select metabolite_id, metabolite_id from hmdb.reaction_elements where reaction_id = %d AND type = 'RightReactionElement'");	
		resultMap.put("GetEnzyme_FromProteinTable_ByReactionID", "select * from hmdb.proteins where id in (select protein_id from hmdb.reaction_enzymes where reaction_id = %d)");
		return resultMap;
	}
}
