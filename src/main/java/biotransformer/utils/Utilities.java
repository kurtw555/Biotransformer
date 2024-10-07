/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.netlib.lapack.Sgbtf2;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.transformation.Biotransformation;
import biotransformer.validateModels.InValidSMARTS;
import biotransformer.validateModels.IsEndProductSMARTS;

public class Utilities {
	public static InValidSMARTS invalidSMARTS = new InValidSMARTS();
	public static IsEndProductSMARTS isEndProductSMARTS = new IsEndProductSMARTS();
	public static SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
	
	public Utilities() {
		// TODO Auto-generated constructor stub
	
	}
	
	public static ArrayList<String> removeDuplicateStrings(ArrayList<String> listOfStrings){
		ArrayList<String> unique =  new ArrayList<String>();		
		LinkedHashSet<String> set = new LinkedHashSet<String>(listOfStrings);
		unique = new ArrayList<String>(set);
		
		return unique;
	}
	
	public static void print(ArrayList<String> aList){
		for(String i : aList){
			System.out.println(i);
		}
	}

	public static IAtomContainerSet getCDKAtomContainerSet(){
		return DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
	}
	
	public static void annotateAtomContainerWithProps(IAtomContainer molecule, LinkedHashMap<Object, Object> props){
		LinkedHashMap<Object, Object> p = new LinkedHashMap<Object, Object>();
		for(Entry<Object, Object> m: props.entrySet()){
			if(m.getKey() != "Synonym"){
				p.put(m.getKey(), m.getValue());
			}
			else {
				p.put(CDKConstants.TITLE,m.getValue());
			}
		}		
		molecule.setProperties(p);	
	}

	/**
	 * Here we select unique biotransformations returning the 
	 * @param biotransformations
	 * @return
	 */
	public static ArrayList<Biotransformation> selectUniqueBiotransformations(ArrayList<Biotransformation> biotransformations){
		ArrayList<Biotransformation> result_bt_list = new ArrayList<>();
		HashMap<String, ArrayList<String>> substrate_metabolite_inChiKey_map = new HashMap<>();
		for(int i = 0; i < biotransformations.size(); i++) {
			boolean duplicate = true;
			Biotransformation bt = biotransformations.get(i);
			//As far as I know, there is no biotransformation contains more than one substrate in this program
			String substrateinChIKey = bt.getSubstrates().getAtomContainer(0).getProperty("InChIKey");
			if(!substrate_metabolite_inChiKey_map.containsKey(substrateinChIKey)) {
				ArrayList<String> metabolite_inChIKey_list = new ArrayList<>();
				for(int j = 0; j < bt.getProducts().getAtomContainerCount(); j++) {
					String oneMetabolite_inChIKey = bt.getProducts().getAtomContainer(j).getProperty("InChIKey");
					if(!metabolite_inChIKey_list.contains(oneMetabolite_inChIKey)) metabolite_inChIKey_list.add(oneMetabolite_inChIKey);
				}
				substrate_metabolite_inChiKey_map.put(substrateinChIKey, metabolite_inChIKey_list);
				duplicate = false;
			}
			else {
				ArrayList<String> metabolite_inChIKey_list = substrate_metabolite_inChiKey_map.get(substrateinChIKey);
				for(int j = 0; j < bt.getProducts().getAtomContainerCount(); j++) {
					String oneMetabolite_inChIKey = bt.getProducts().getAtomContainer(j).getProperty("InChIKey");
					if(!metabolite_inChIKey_list.contains(oneMetabolite_inChIKey)) {
						metabolite_inChIKey_list.add(oneMetabolite_inChIKey);
						duplicate = false;
					}									
				}
				substrate_metabolite_inChiKey_map.put(substrateinChIKey, metabolite_inChIKey_list);	
			}
			if(!duplicate) result_bt_list.add(bt);
		}
		return result_bt_list;
	}
	
	public static ArrayList<Biotransformation> selectUniqueBiotransformations_old(ArrayList<Biotransformation> biotransformations){
		ArrayList<Biotransformation> unique_bts = new ArrayList<Biotransformation>();
		for(int i = 0; i < biotransformations.size(); i ++){
			if(!containsBiotransformation(unique_bts, biotransformations.get(i))){
				unique_bts.add(biotransformations.get(i));
			}
		}	
		try {
			unique_bts = selectByUniqueMetabolites(unique_bts);
		}catch (Exception e) {
			return unique_bts;
		}
		return unique_bts;
	}
	
	public static boolean containsBiotransformation(ArrayList<Biotransformation> biotransformations, Biotransformation bt){
		boolean inc = false;
		for(int i = 0; i < biotransformations.size(); i ++){
			if(biotransformations.get(i).equals(bt)){
				inc = true;
				break;
			}
		}
		return inc;
	}
	public static IAtomContainerSet createEmptyAtomContainerSet(){
		return DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	}
	
	
	public static String returnFirstCleanSynonym(String[] synonyms){
		String fcs = null;
		Pattern p = Pattern.compile("CHEBI:[0-9]+|UNII-|CHEMBL|ZINC|DB[0-9]+|[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-|[0-9]+-[0-9]+-[0-9]|^AC[0-9]+|%",  Pattern.CASE_INSENSITIVE);
		for(int i=0; i < synonyms.length ; i++){
			Matcher b = p.matcher(synonyms[i]);
			if(!b.find()){
				fcs = synonyms[i];
				break;				
			}		
		}
		if(fcs == null){
			fcs = synonyms[0];
		}
		return fcs;		
	}

	public static String returnFirstCleanSynonym(ArrayList<String> synonyms){
		String fcs = null;
		// CHEMSPIDER|CID[0-9]+|SID[0-9]+|
		Pattern p = Pattern.compile("CHEBI:[0-9]+|UNII-|CHEMBL|ZINC|DB[0-9]+|[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-|[0-9]+-[0-9]+-[0-9]|^AC[0-9]+|%",  Pattern.CASE_INSENSITIVE);
		for(int i=0; i < synonyms.size() ; i++){
			Matcher b = p.matcher(synonyms.get(i).trim());			
			if(!b.find()){
				fcs = synonyms.get(i);
				break;			
			}		
		}
		if(fcs == null){
			fcs = synonyms.get(0);
		}
		return fcs;	
	}

	public static void addPhysicoChemicalProperties(IAtomContainer molecule) throws CDKException {			
			LinkedHashMap<String, String> properties = ChemStructureExplorer.computePhysicoChemicalProperties(molecule);			
			for(Map.Entry<String, String> prop : properties.entrySet()){
				molecule.setProperty(prop.getKey(), String.format("%.8s", Double.valueOf(prop.getValue()))   );
			}
			
		}

	/**
	 * Added by Siyang. To select unique biotransformations by metabolites
	 * @param biotransformations
	 * @return
	 */
	public static ArrayList<Biotransformation> selectByUniqueMetabolites(ArrayList<Biotransformation> biotransformations) throws Exception{
		ArrayList<Biotransformation> unique_bts = new ArrayList<>();
		ArrayList<String> stored_InChIKey = new ArrayList<>();
		for(int i = 0; i < biotransformations.size(); i++){
			boolean isDuplicate= true;
			IAtomContainerSet metabolites = biotransformations.get(i).getProducts();
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				IAtomContainer oneMetabolite = metabolites.getAtomContainer(j);
				String inchikey = oneMetabolite.getProperty("InChIKey");
				if(!stored_InChIKey.contains(inchikey)){
					stored_InChIKey.add(inchikey);
					isDuplicate = false;
				}
			}
			if(isDuplicate) continue;
			else unique_bts.add(biotransformations.get(i));
		}
		return unique_bts;
	}

	public static IAtomContainer create3DCoordinates(IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(oneMole.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(oneMole);
		adder.addImplicitHydrogens(oneMole);
		IChemObjectBuilder builder_1 = SilentChemObjectBuilder.getInstance();
		ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder_1);
		IAtomContainer oneMolecule = mb3d.generate3DCoordinates(oneMole, false);
		AtomContainerManipulator.suppressHydrogens(oneMolecule);
		return oneMolecule;
	}
	public static IAtomContainerSet extractProductsFromBiotransformations(ArrayList<Biotransformation> biotransformations) throws Exception{
		IAtomContainerSet acontainers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		LinkedHashMap<String, IAtomContainer> hMap = new LinkedHashMap<String, IAtomContainer>();
		for(Biotransformation b : biotransformations){
			int numAtoms_substrate = 0;
			IAtomContainerSet substrates = b.getSubstrates();
			for(int i = 0; i < substrates.getAtomContainerCount(); i++){
				IAtomContainer oneSubstrate = substrates.getAtomContainer(i);
				//AtomContainerManipulator.suppressHydrogens(oneSubstrate);
				//numAtoms_substrate += oneSubstrate.getAtomCount(); 
				numAtoms_substrate += countNonHydrogenAtoms(oneSubstrate);
			}
			for(IAtomContainer ac : b.getProducts().atomContainers()){
				String ikey = ac.getProperty("InChIKey");			
				//The suppressHydrogens function may not work for ac. The assumption here is not that useful and doesn't make sense
				int numAtoms_metabolites = countNonHydrogenAtoms(ac);
				//Use 5 here considering the number of non-hydrogen atoms in sulfate (was glycine)
				if(numAtoms_metabolites >= numAtoms_substrate + 4 || ChemStructureExplorer.getMajorIsotopeMass(ac) > 1500){
					ac.setProperty("isEndProduct", true);
				}
				if(!hMap.containsKey(ikey)){
					hMap.put(ikey, ac);
					
				}
			}
		}
		ArrayList<IAtomContainer> at =  new ArrayList<IAtomContainer>( hMap.values());
//		System.err.println(at.get(0).getClass());
		for(IAtomContainer a : at){
			acontainers.addAtomContainer(ChemStructureManipulator.preprocessContainer(a));
		}
		return acontainers;
//		return ChemStructureExplorer.uniquefy(acontainers);
	}
	/**
	 * This 
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public static int countNonHydrogenAtoms(IAtomContainer molecule) throws Exception{
		int count = 0;
		for(int i = 0; i < molecule.getAtomCount(); i++){
			if(!molecule.getAtom(i).getSymbol().equalsIgnoreCase("H")) count++;
		}
		return count;
	}
	/**
	 * Here we don't want to consider the metabolites that have too few atoms in the structure.
	 * @param oneMole
	 * @return
	 */
	public static boolean isValidMetabolte(IAtomContainer oneMole) throws Exception{
		
		if(invalidSMARTS.containInvalidSubstructure(oneMole)) return false;
		int countAtom_nonH = 0;
		int countCarbon = 0;
		for(int t = 0; t < oneMole.getAtomCount(); t++){
			if(!oneMole.getAtom(t).getSymbol().equalsIgnoreCase("H")) countAtom_nonH++;
			if(oneMole.getAtom(t).getSymbol().equalsIgnoreCase("C")) countCarbon++;
		}
		if(countAtom_nonH >=1 && countCarbon >= 1) return true;
		else return false;
	}
	
	public static ArrayList<String> updateProcessedSubstratePool(ArrayList<String> substrateProcessed, IAtomContainerSet novelSubstrateSet) throws Exception{
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();			
		for(int j = 0; j < novelSubstrateSet.getAtomContainerCount(); j++){
			IAtomContainer oneSubstrate = novelSubstrateSet.getAtomContainer(j);
			String subStrate_inChiKey = inchiGenFactory.getInChIGenerator(oneSubstrate).getInchiKey();
			if(!substrateProcessed.contains(subStrate_inChiKey)) substrateProcessed.add(subStrate_inChiKey);
		}
		return substrateProcessed;
	}
	public static ArrayList<String> updateProcessedSubstratePool(ArrayList<String> substrateProcessed, IAtomContainer novelSubstrate) throws Exception{
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();					
		IAtomContainer oneSubstrate = novelSubstrate;
		String subStrate_inChiKey = inchiGenFactory.getInChIGenerator(oneSubstrate).getInchiKey();
		if(!substrateProcessed.contains(subStrate_inChiKey)) substrateProcessed.add(subStrate_inChiKey);
		
		return substrateProcessed;
	}
	/**
	 * This function will explore the substratePool and return the substrates that have not been input into the metabolism iteration
	 * @param substrateSet
	 * @param processedInChiKeyList
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet getNovelSubstrates(IAtomContainerSet substrateSet, ArrayList<String> processedInChiKeyList) throws Exception{
		ArrayList<String> uniqueInChIKeyList= new ArrayList<>();
		IAtomContainerSet cleanedSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		for(int i = 0; i < substrateSet.getAtomContainerCount(); i++){
			IAtomContainer oneMole = substrateSet.getAtomContainer(i);
			//if(this.isEndProductSMARTS.containEndSubstructure(oneMole)) continue;
			if(isEndProductSMARTS.isEndProduct(oneMole)) continue;			
			InChIGenerator inchiGen = inchiFactory.getInChIGenerator(oneMole);
			String inChiKey = inchiGen.getInchiKey();
			if(!processedInChiKeyList.contains(inChiKey) && !uniqueInChIKeyList.contains(inChiKey)){
				uniqueInChIKeyList.add(inChiKey);
				cleanedSet.addAtomContainer(oneMole);
			}
		}
		return cleanedSet;
	}
}
