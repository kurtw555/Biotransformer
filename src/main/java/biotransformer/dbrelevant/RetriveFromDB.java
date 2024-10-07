package biotransformer.dbrelevant;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.Utilities;
import biotransformer.utils.filterCertainClasses.FilterCertainClasses;

/**
 * Because most reactions/enzymes in HMDB are EC based and the PhaseI/PhaseII/GutMicrobial information are not given, this RetriveFromHMDB class is used as a module in ECBased biotransformation prediction
 * query_biotransformation_type = all, cyp450, phaseII, microbial. Please note that all = ecb in this case.
 * @author Tian
 *
 */
public class RetriveFromDB {
	private final String DB_URL;
	private final String USER;
	private final String PASS;
	private ArrayList<String> inChIKeyList_ignore;
	private boolean first14 = true;
	private boolean checkOverflow = true;
	private boolean substituted = false;
	private boolean useSubstituition = false;
	private int overflow_threshold;
	public String data_source;
	public FilterCertainClasses fcc = new FilterCertainClasses();
	RetrieveFromDBHelpingFunctions rfdb_help;
	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	//private String query_biotransformation_type;
	public static void main(String[] args) throws Exception{
		//String source = "pathbank";
		String source = "hmdb";		
		boolean checkExported = true;
		boolean checkOverflow = true;
		int overflow_threshold = 30;
		String query_biotransformation_type = "all";
		RetriveFromDB rfdb = new RetriveFromDB(source, query_biotransformation_type, checkOverflow, overflow_threshold, true);
		//String InChIKey = "MUMGGOZAMZWBJJ-UHFFFAOYSA-N";
		//String InChIKey = "PFTAWBLQPZVEMU-UKRRQHHQSA-N";
		//String InChIKey = "WQZGKKKJIJFFOK-GASJEMHNSA-N";// ? D-glucose?
		
		//String InChIKey = "WQZGKKKJIJFFOK-VFUOTHLCSA-N"; //D-glucose
		//String InChIKey = "NBSCHQHZLSJFNQ-VFUOTHLCSA-N";//Glucose-6-phosphate
		//String InChIKey = "BGWGXPAPYGQALX-ARQDHWQXSA-N";//Fructose-6-phosphate
		//String InChIKey = "RNBGYGVWRKECFJ-ARQDHWQXSA-N";//Fructose-1,6-biphosphate
		//String InChIKey = "LXJXRIRHZLFYRP-VKHMYHEASA-N";//G3P
		//String InChIKey = "LJQLQCAXBUHEAZ-UWTATZPHSA-N";//1,3-biphosphoglycerate
		//String InChIKey = "OSJPPGNTCRNQQC-UHFFFAOYSA-N";//3-phosphoglycerate
		//String InChIKey = "GXIURPTVHJPJLF-UWTATZPHSA-N";//2-phosphoglycerate
		String InChIKey = "DTBNBXWJWCWCIK-UHFFFAOYSA-N"; //Phosphoenolpyruvate
		
		//String InChIKey = "MNBKLUUYKPBKDU-BBECNAHFSA-N"; //lipid
		//String InChIKey = "QIVBCDIJIAJPQS-VIFPVBQESA-N";//Amino acid
		//String InChIKey = "LCTONWCANYUPML-UHFFFAOYSA-N";//Pyruvate
		//String InChIKey = "MTCFGRXMJLQNBG-REOHCLBHSA-N";//Serine
		//String InChIKey = "RGJOEKWQDUBAIZ-IBOSZNHHSA-N";
		//String InChIKey = "FBPFZTCFMRRESA-GUCUJZIJSA-Na";
		//String InChIKey = "WQZGKKKJIJFFOK-ZZWDRFIYSA-N";//L-glucose
		//String InChIKey = "GHOKWGTUZJEAQD-ZETCQYMHSA-N";//Not-Endogenous
		//ArrayList<Biotransformation> result = rfdb.getPredictedBiotransformation(InChIKey, checkExported);
		//System.out.println(rfdb.isEndogenous(InChIKey));
		
		Biotransformer bf = new Biotransformer(BioSystemName.HUMAN, false, false);
		ArrayList<Biotransformation> result = rfdb.getPredictedBiotransformationContainsElement("Se"); 
		bf.saveBioTransformationProductsToSdf(result, "G:/Biology/SanityTest/Selenium.sdf");
		//ArrayList<Biotransformation> result = rfdb.getPredictedBiotransformationContainsElement("Sulfate");
		//bf.saveBioTransformationProductsToSdf(result, "G:/Biology/SanityTest/Sulfate.sdf");
		//ArrayList<Biotransformation> result = rfdb.getPredictedBiotransformationContainsElement("Phosphate");
		//bf.saveBioTransformationProductsToSdf(result, "G:/Biology/SanityTest/Phosphate.sdf");
		System.out.println(result.size());
	}
	
	public RetriveFromDB(String source, String query_biotransformation_type, boolean checkOverflow, int overflow_threshold, boolean useSubstituition) throws Exception{
		//this.query_biotransformation_type = query_biotransformation_type;
		this.useSubstituition = useSubstituition;
		if(source.equalsIgnoreCase("pathbank")){
			this.DB_URL = "jdbc:mysql://192.241.161.227/pathbank";
			this.USER = "";
			this.PASS = "";
		}
		else{			
			this.DB_URL = "jdbc:mysql://162.243.63.130:3306/hmdb";
			this.USER = "biotransformer";
			this.PASS = "cyp450";
		}
		this.checkOverflow = checkOverflow;
		this.overflow_threshold = overflow_threshold;
		this.data_source = source;
		this.rfdb_help = new RetrieveFromDBHelpingFunctions(source, query_biotransformation_type, this.first14, this.checkOverflow, this.substituted, this.overflow_threshold);
		this.inChIKeyList_ignore = rfdb_help.getInChIKeyList_ignore();
		
	}
	public RetriveFromDB(String query_biotransformation_type, boolean useSubstituition) throws Exception{
		//this.query_biotransformation_type = query_biotransformation_type;
		this.useSubstituition = useSubstituition;
		this.DB_URL = "jdbc:mysql://162.243.63.130:3306/hmdb";
		this.USER = "biotransformer";
		this.PASS = "cyp450";
		this.checkOverflow = true;
		this.overflow_threshold = 30;//by default
		this.rfdb_help = new RetrieveFromDBHelpingFunctions("hmdb", query_biotransformation_type,  this.first14, this.checkOverflow, this.substituted, this.overflow_threshold);
		this.inChIKeyList_ignore = rfdb_help.getInChIKeyList_ignore();
		this.data_source = "hmdb";//By default
		
	}
	/**
	 * This function allows will retrieve metabolites from pathbank for the input set of molecules for a certain number of steps.
	 * Note that if the compound is retrieved from pathBank, then itself and all its future metabolites will be retrieved from Pathbank
	 * @param molecules
	 * @param steps
	 * @param checkExported
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getBiotransformationRetrievedFromDB_setAndstep(IAtomContainer molecule, boolean checkExported) throws Exception{
		ArrayList<Biotransformation> resultArray = new ArrayList<>();
		InChIGeneratorFactory inChIGenerator = InChIGeneratorFactory.getInstance();
		InChIGenerator inchiGen; 
		String inChiKey;			
		ArrayList<Biotransformation> tempResultArray = new ArrayList<>();			
		IAtomContainer oneMole = molecule;
		inchiGen = inChIGenerator.getInChIGenerator(oneMole);
		inChiKey = inchiGen.getInchiKey();
		ArrayList<Biotransformation> oneResultArray = getBiotransformationRetrievedFromDB(oneMole, checkExported);
		tempResultArray.addAll(oneResultArray);							
		//IAtomContainerSet tempProducts = this.ecb.extractProductsFromBiotransformations(tempResultArray);			
		resultArray.addAll(tempResultArray);
		
		System.out.println("Finished DB process");
		return Utilities.selectUniqueBiotransformations(resultArray);
	}
	/**
	 * This function allows will retrieve metabolites from pathbank for the input set of molecules for a certain number of steps.
	 * Note that if the compound is retrieved from pathBank, then itself and all its future metabolites will be retrieved from Pathbank
	 * @param molecules
	 * @param steps
	 * @param checkExported
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getBiotransformationRetrievedFromDB_setAndstep_set(IAtomContainerSet molecules, int steps, boolean checkExported) throws Exception{
		ArrayList<Biotransformation> resultArray = new ArrayList<>();
		ArrayList<String> processedInChIKeyArray = new ArrayList<>();
		InChIGeneratorFactory inChIGenerator = InChIGeneratorFactory.getInstance();
		InChIGenerator inchiGen; 
		String inChiKey;
		IAtomContainerSet substarteSet = molecules;
		while(steps > 0){
			steps = steps - 1;			
			ArrayList<Biotransformation> tempResultArray = new ArrayList<>();
			for(int i = 0; i < substarteSet.getAtomContainerCount(); i++){
				IAtomContainer oneMole = substarteSet.getAtomContainer(i);
				inchiGen = inChIGenerator.getInChIGenerator(oneMole);
				inChiKey = inchiGen.getInchiKey();
				if(processedInChIKeyArray.contains(inChiKey)) continue;
				else processedInChIKeyArray.add(inChiKey);
				ArrayList<Biotransformation> oneResultArray = getBiotransformationRetrievedFromDB(oneMole, checkExported);
				//If no metabolites retrieved from HMDB, then we use Biotransformer to make predictions for the query molecule.
				if(oneResultArray == null || oneResultArray.isEmpty()) {
					oneResultArray = predictBiotransformation(oneMole);
				}
				tempResultArray.addAll(oneResultArray);				
			}
			IAtomContainerSet tempProducts = Utilities.extractProductsFromBiotransformations(tempResultArray);
			substarteSet = tempProducts;
			
			resultArray.addAll(tempResultArray);
		}
		System.out.println("Finished DB process");
		return Utilities.selectUniqueBiotransformations(resultArray);
	}
	/**
	 * This is the function that should be called in the ECBased metabolism prediction
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getBiotransformationRetrievedFromDB(IAtomContainer molecule, boolean checkExported) throws Exception{
		InChIGeneratorFactory inChIGenerator = InChIGeneratorFactory.getInstance();
		InChIGenerator inchiGen = inChIGenerator.getInChIGenerator(molecule);
		String inChIKey = inchiGen.getInchiKey();
		String first14 = inChIKey.split("-")[0];
		//Return empty arraylist if the query molecule's inChiKey's first 14 characters are in the ignore list.
		if(this.inChIKeyList_ignore.contains(first14)) return new ArrayList<Biotransformation>();
		if(inChIKey.equals("WQZGKKKJIJFFOK-GASJEMHNSA-N") || inChIKey.equals("WQZGKKKJIJFFOK-DVKNGEFBSA-N") || inChIKey.equals("WQZGKKKJIJFFOK-UHFFFAOYNA-N") || inChIKey.equals("WQZGKKKJIJFFOK-UHFFFAOYSA-N")) inChIKey = "WQZGKKKJIJFFOK-VFUOTHLCSA-N";
		ArrayList<Biotransformation> results = getPredictedBiotransformation(inChIKey, checkExported);
		return results;
	}
	/**
	 * This function will be used in the ruby rails app and return ArrayList<Biotransformation>
	 * @param inChIKey
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getPredictedBiotransformation(String inChiKey, boolean checkExported) throws Exception{
		String first14 = inChiKey.split("-")[0];
		//Return empty arraylist if the query molecule's inChiKey's first 14 characters are in the ignore list.
		if(this.inChIKeyList_ignore.contains(first14)) return new ArrayList<Biotransformation>();
		ArrayList<Biotransformation> resultArrayList = new ArrayList<>();
		try(
				Connection conn = DriverManager.getConnection(this.DB_URL, this.USER, this.PASS);
				Statement stmt = conn.createStatement();
			    ){		      
					ArrayList<Biotransformation> oneResult = this.rfdb_help.getBiotransformationInfo(stmt, inChiKey, checkExported, this.useSubstituition);
					resultArrayList.addAll(oneResult);
			      } catch (Exception e) {
					e.printStackTrace();
				}
		//System.out.println(resultArrayList.size());
		return resultArrayList;
	}

	
	/**
	 * This function checks if the molecule is an endogenous or not.
	 * If it's an endogenous molecule, then retrieve the metabolite from database.
	 * Otherwise, use Biotransformer to make predictions.
	 * This function is used to replace the retriveFromPathbank function which means we don't classify it by the molecule classes but if it's an endogeous instead.
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public boolean isEndogenous(String inChIKey) throws Exception{
		String query = "select metabolites.id, metabolites.moldb_inchikey, ontologies.origin from metabolites, ontologies where ontologies.metabolite_id = metabolites.id and export_to_hmdb =1 and moldb_inchikey like '%s%%'"; //origin='Endogenous'
		
		if(this.useSubstituition) {
			query = String.format(query, inChIKey.split("-")[0]);
		}
		else{//We now check the full InChIkey instead of the first 14 characters to determine if the molecule is endogenous or not
			query = String.format(query, inChIKey);
		}
		//"select id, export_to_hmdb, moldb_inchikey, moldb_inchi, moldb_smiles, moldb_average_mass, moldb_mono_mass, moldb_formula, moldb_alogps_logp as moldb_alogps_logs from hmdb.metabolites where moldb_inchikey = '%s'";
		try(
				Connection conn = DriverManager.getConnection(this.DB_URL, this.USER, this.PASS);
				Statement stmt = conn.createStatement();
			    ){		      
					ResultSet rs = stmt.executeQuery(query);
					while(rs.next()){					
						String origin = rs.getString("origin");
						if(origin==null) continue;
						if(origin.equals("Endogenous") || origin.equals("endogenous")){
							//System.out.println(rs.getString("moldb_inchikey"));
							return true;							
						}
						if(origin.contains("Endogenous") || origin.contains("endogenous")) {
							if(origin.contains("Food") || origin.contains("food")) {
								//System.out.println(rs.getString("moldb_inchikey"));
								return true;	
							}
						}
					}
					return false;
					//System.out.println(origin);
			      } catch (Exception e) {
					e.printStackTrace();
				}
		return false;	
	}
	
	
	/**
	 * This function will be used in the ruby rails app and return ArrayList<Biotransformation>
	 * @param inChIKey
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getPredictedBiotransformationContainsElement(String keyElement) throws Exception{
		ArrayList<Biotransformation> resultArrayList = new ArrayList<>();
		HashMap<String, Object> precursorMap = getPrecursorMap("Se");
		String enzyme = "Unkown";
		String reaction = "Unkown";
		String queryKeyWord = "Se";
		if(keyElement.equals("Sulfate")) {
			enzyme = "SULTs";
			reaction = "Sulfation";
			queryKeyWord = "S";
		}
		else if(keyElement.equals("Phosphate")) {
			enzyme = "PTS";
			reaction = "Phosphation";
			queryKeyWord = "%P";
		}
		try(
				Connection conn = DriverManager.getConnection(this.DB_URL, this.USER, this.PASS);
				Statement stmt = conn.createStatement();
			    ){		      	
					String query = this.rfdb_help.queryTable.get("GetMetabolitesContainsElements");
					query = String.format(query, "%" + queryKeyWord);
					ResultSet rs = stmt.executeQuery(query);
					HashMap<String, Object> oneMetaboliteMap = new HashMap<>();
					while(rs.next()) {
						String moldb_smiles = rs.getString("moldb_smiles");
						if(!structContainsFunc(moldb_smiles, keyElement)) continue;
						Double moldb_average_mass = rs.getDouble("moldb_average_mass");
						Double moldb_mono_mass = rs.getDouble("moldb_mono_mass");
						String moldb_formula = rs.getString("moldb_formula");
						Double moldb_alogps_logp = rs.getDouble("moldb_alogps_logs");
						String moldb_inchikey = rs.getString("moldb_inchikey");
						String moldb_inchi = rs.getString("moldb_inchi");
						//Close this result set every time after it is used			
						oneMetaboliteMap.put("moldb_smiles", moldb_smiles);
						oneMetaboliteMap.put("moldb_average_mass", moldb_average_mass);
						oneMetaboliteMap.put("moldb_mono_mass", moldb_mono_mass);
						oneMetaboliteMap.put("moldb_formula", moldb_formula);
						oneMetaboliteMap.put("moldb_alogps_logp", moldb_alogps_logp);
						oneMetaboliteMap.put("moldb_inchikey", moldb_inchikey);
						oneMetaboliteMap.put("moldb_inchi", moldb_inchi);											
						//System.out.println(moldb_smiles);
						if(this.rfdb_help.generationBiotransformation(precursorMap, oneMetaboliteMap, enzyme, reaction) != null) resultArrayList.add(this.rfdb_help.generationBiotransformation(precursorMap, oneMetaboliteMap, enzyme, reaction));
						
					}
					
			      } catch (Exception e) {
					e.printStackTrace();
				}
		//System.out.println(resultArrayList.size());
		return resultArrayList;
	}
	
	
	/**
	 * This function generates the biotranformer predicted results for those "HMDB valid" compound that doesn't have any returned metabolites retrieved from HMDb 
	 * all, cyp450, phaseII, microbial. Note all = ecb
	 * @param substrate
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> predictBiotransformation(IAtomContainer substrate) throws Exception{
		return new ArrayList<>();
//		if(this.query_biotransformation_type.equals("all")) {
//			return this.hsbt.simulateHumanMetabolism(substrate, 1);
//		}
//		else if(this.query_biotransformation_type.equals("cyp450")) {
//			return this.cyp450bt.predictCyp450BiotransformationChainByMode(substrate, true, true, 1, 0.0, 3);
//		}
//		else if(this.query_biotransformation_type.equals("phaseII")) {
//			return this.phase2b.applyPhase2TransformationsChainAndReturnBiotransformations(substrate, true, true, true, 1, 0.0);
//		}
//		else if(this.query_biotransformation_type.equals("microbial")) {
//			return this.hgut.simulateGutMicrobialMetabolism(substrate, true, true, 1, 0.0);
//		}
//		else throw new Exception("Please specify the biotransformation type you want to query when calling retrieving from db feature. The types are all, cyp450, phaseII and microbial");
	}
	
	
	public HashMap<String, Object> getPrecursorMap(String type) {
		HashMap<String, Object> resultMap = new HashMap<>();
		String moldb_smiles = "[Se++]";
		Double moldb_average_mass = 78.960;
		Double moldb_mono_mass = 79.917;
		String moldb_formula = "Se";
		Double moldb_alogps_logp = 5.0;
		String moldb_inchikey = "MFSBVGSNNPNWMD-UHFFFAOYSA-N";
		String moldb_inchi = "InChI=1S/Se/q+2";	
		if(type.contains("Se")) {
			
		}
		else if(type.contains("Sulfate")) {
			moldb_smiles = "OS(O)(=O)=O";
			moldb_average_mass = 98.07;
			moldb_mono_mass = 97.967;
			moldb_formula = "H2O4S";
			moldb_alogps_logp = -1.4;
			moldb_inchikey = "QAOWNCQODCNURD-UHFFFAOYSA-N";
			moldb_inchi = "InChI=1S/H2O4S/c1-5(2,3)4/h(H2,1,2,3,4)";	
		}
		else if(type.contains("Phosphate")) {
			moldb_smiles = "OP(O)(O)=O";
			moldb_average_mass = 97.99;
			moldb_mono_mass = 97.98;
			moldb_formula = "H3O4P";
			moldb_alogps_logp = -2.1;
			moldb_inchikey = "NBIIXXVUZAFLBC-UHFFFAOYSA-N";
			moldb_inchi = "InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)";	
		}
		resultMap.put("moldb_smiles", moldb_smiles);
		resultMap.put("moldb_average_mass", moldb_average_mass);
		resultMap.put("moldb_mono_mass", moldb_mono_mass);
		resultMap.put("moldb_formula", moldb_formula);
		resultMap.put("moldb_alogps_logp", moldb_alogps_logp);
		resultMap.put("moldb_inchikey", moldb_inchikey);
		resultMap.put("moldb_inchi", moldb_inchi);	
		
		return resultMap;
	}
	
	public boolean structContainsFunc(String smiles, String type) throws Exception{
		IAtomContainer molecule = this.sp.parseSmiles(smiles);
		String smarts = "[#8]S([#8])(=O)=O";
		if(type.equals("Se")) return true;
		if(type.equals("Phosphate")) {
			smarts = "[#8]P([#8])([#8])=O";
		}
		SmartsPattern smp = SmartsPattern.create(smarts);
		if(smp.matches(molecule)) return true;
		else return false;
	}
}
