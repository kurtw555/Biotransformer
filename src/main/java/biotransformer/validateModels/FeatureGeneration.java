package biotransformer.validateModels;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;


public class FeatureGeneration {
	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	/**
	 * definition for all the adopted SYBYL atom types chosen by us
	 * 23 in total
	 */
	public static final String[] atomTypeLookupTable = { "C.1", "C.2", "C.3", "C.ar", "C.cat", "N.1",
			"N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2",
			"S.3", "S.O", "S.O2", "P.3", "F", "Cl", "Br", "I" };
	public static final String[] featuresToRemove_Urine = {"SMILES", "C.cat" , "N.1", "F", "Cl", "Br", "I", "halogens", "nitriles", "dihydropyridines", "alkyl carbamates", "sulfonamides" , "hydrazone groups"};
	public static final String[] featuresToRemove_ECBased = {"SMILES", "C.cat", "N.1", "F", "Cl", "I", "nitriles", "dihydropyridines", "hydrazone groups"};
	public static final String[] featuresToRemove_gut_substrate = {"SMILES", "C.1", "C.cat", "N.1", "F", "Cl", "Br", "I", "halogens", "nitriles", "alkyl carbamates", "sulfone groups", "tert-alicyclic amines", "terminal acetylenes", "morpholine rings", "sulfonamides", "hydrazone groups", "thiophene rings"};
	public static final String[] featuresToRemove_gut_metabolite = {"SMILES", "C.cat", "S.O", "F", "Cl", "Br", "I", "halogens", "alkyl carbamates", "sulfone groups", "tert-alicyclic amines", "terminal acetylenes", "morpholine rings", "sulfonamides", "thiophene rings", "hydrazone groups"};
	public static final String[] continuousFeatures = {"mass", "rotatableBonds",	"acceptor",	"donor", "logD", "tpsa"};
	//SMILES string is an identifier and should not be used as a feature
	

	/**
	 * This function will generate 
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public String generateMolecularFeature(IAtomContainer oneMole) throws Exception{
		String molecularFeature = generateMolecularPropertyFeature(oneMole);
		String molecularAtomTypeFeature = getAtomTypeFeatures(oneMole);
		String fingerprintFeature = generateFingerPrintFeature(oneMole);
		String resultString = molecularFeature + "," + molecularAtomTypeFeature + "," + fingerprintFeature;
		return resultString;
		//molecularFeature.split(",");
		//fingerprintFeature.split(",");
	}
	
	/**
	 * This function will generate the title line for the data file that contains all feature values
	 * Please note that we add the SMILES as an identifier in the feature value data file.
	 * @return
	 */
	public String generateTitleLine(){
		String molecularProperty = "mass,rotatableBonds,acceptor,donor,logD,tpsa";
		StringBuffer sb_atomType = new StringBuffer();
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(i == atomTypeLookupTable.length-1){
				sb_atomType.append(atomTypeLookupTable[i]);
			}
			else{
				sb_atomType.append(atomTypeLookupTable[i]).append(",");
			}
		}
		String atomTypeFeature_title = sb_atomType.toString();
		StringBuffer sb_fingerPrint = new StringBuffer();
		HashMap<String, String> smarts_pattern_map = getSMARTS_Pattern();
		for(String name : smarts_pattern_map.keySet()){
			sb_fingerPrint.append(name).append(",");
		}
		sb_fingerPrint.deleteCharAt(sb_fingerPrint.length()-1);
		String atomFingerPrint_title = sb_fingerPrint.toString();
		
		return "SMILES" + "," + molecularProperty + "," + atomTypeFeature_title + "," + atomFingerPrint_title;
	}
	/**
	 * This function will generate molecular features for the query molecule
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public String generateMolecularPropertyFeature(IAtomContainer oneMole) throws Exception{		
		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(oneMole);		
		double mass = MolecularFormulaManipulator.getMajorIsotopeMass(formula);
		
		RotatableBondsCountDescriptor Ro = new RotatableBondsCountDescriptor();
		double rotatableBonds = Double.parseDouble(Ro.calculate(oneMole).getValue().toString());
		
		HBondAcceptorCountDescriptor acceptor_des = new HBondAcceptorCountDescriptor();
		double acceptor = Double.parseDouble(acceptor_des.calculate(oneMole).getValue().toString());
		
		HBondDonorCountDescriptor donor_bonds = new HBondDonorCountDescriptor();
		double donor = Double.parseDouble(donor_bonds.calculate(oneMole).getValue().toString());
		
		//PredictLogD logD_predictor = new PredictLogD();
		//Double logD = logD_predictor.predictLogD(oneMole);
		
		TPSADescriptor tpsa_des = new TPSADescriptor();
		String tpsa_value = tpsa_des.calculate(oneMole).getValue().toString();
		Double tpsa = Double.parseDouble(tpsa_value);
		
		String resultString = mass + "," + rotatableBonds + "," + acceptor + "," + donor + "," + 0.0 + "," + tpsa;
		return resultString;
	}
	/**
	 * This function will generate the molecular fingerprint features for the query molecule
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public String generateFingerPrintFeature(IAtomContainer oneMole) throws Exception{
		HashMap<String, String> smartsPattern_map = getSMARTS_Pattern();
		//IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		StringBuffer sb = new StringBuffer();
		for(String name : smartsPattern_map.keySet()){
			String smirts = smartsPattern_map.get(name);
			//Pattern pattern = SmartsPattern.create(smirts);
			Pattern pattern = SmartsPattern.create(smirts);
			int occurrence = pattern.matchAll(oneMole).count();
			sb.append(occurrence).append(",");
			//System.out.println(pattern.matchAll(oneMole).count());
			
		}
		sb.deleteCharAt(sb.length()-1);
		return sb.toString();
	}
	
	/**
	 * This function will generate the atom types features the numClosestAtoms for the oneAtom within the molecule
	 * @param oneMole
	 * @param oneAtom
	 * @param numClosestAtoms
	 * @return
	 * @throws Exception
	 */
	public String getAtomTypeFeatures(IAtomContainer molecule) throws Exception{
		ArrayList<IAtom> closestAtoms = new ArrayList<>();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		for(int i = 0; i < molecule.getAtomCount(); i++){
			closestAtoms.add(molecule.getAtom(i));
		}
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		String[] resultStrArray = new String[atomTypeLookupTable.length];	
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}		
		for(int k = 0; k < closestAtoms.size(); k++){
			Integer atomIdx = molecule.indexOf(closestAtoms.get(k));
			IAtomType[] atomTypeList = typeMatcher.findMatchingAtomTypes(molecule); 
			IAtomType atomType = atomTypeList[atomIdx];
			String type;
			if(atomType != null){
				type = atomType.getAtomTypeName();
			}
			else{
				type = molecule.getAtom(atomIdx).getSymbol();
			}
			//Find the matched atomType and update the value to current + "1"
			for(int i = 0; i < atomTypeLookupTable.length; i++){
				if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
					resultStrArray[i] = String.valueOf(Integer.parseInt(resultStrArray[i]) + 1);// "1";
				}
			}		
		}	
		StringBuffer resultStrBuffer = new StringBuffer();
		for(int k = 0; k <resultStrArray.length; k++){
			if(k==0){
				resultStrBuffer.append(resultStrArray[0]);
			}
			else{
				resultStrBuffer = resultStrBuffer.append(",").append(resultStrArray[k]);
			}
		}
		return resultStrBuffer.toString();
	}
	
	
	/**
	 * This function contains the list of SMARTS strings for important functional groups mentioned in the paper:
	 * Mapping human microbiome drug metabolism by gut bacteria and their genes; see extended data Fig. 2
	 * The SMARTS strings are obtained from RDK toolkit's fragmentDescriptor class which can be found at https://github.com/rdkit/rdkit/blob/b208da471f8edc88e07c77ed7d7868649ac75100/Data/FragmentDescriptors.csv#L58 
	 * @return
	 */
	public static HashMap<String, String> getSMARTS_Pattern(){
		HashMap<String, String> resultMap = new HashMap<>();
		resultMap.put("Number of aliphatic hydroxyl groups", "[C!$(C=O)]-[OH]");
		resultMap.put("cyclic esters", "[C&R1](=O)[O&R1][C&R1]");
		resultMap.put("thiazole rings", "c1scnc1");
		resultMap.put("nitriles", "[NX1]#[CX2]");
		resultMap.put("terminal acetylenes", "C#[CH]");
		resultMap.put("furan rings", "o1cccc1");
		resultMap.put("alkyl carbamates", "C[NH1]C(=O)OC");
		resultMap.put("allylic oxidation sites with exclusion", "[$(C=C-C);!$(C=C-C-[N,O,S]);!$(C=C-C-C-[N,O]);!$(C12=CC(=O)CCC1C3C(C4C(CCC4)CC3)CC2)]");
		resultMap.put("carboxylic acids", "[#6]C(=O)[O;H,-1]");
		resultMap.put("carboxylic acids with O2", "[CX3](=O)[OX1H0-,OX2H1]");
		resultMap.put("esters", "[#6][CX3](=O)[OX2H0][#6]");
		resultMap.put("ketones", "[#6][CX3](=O)[#6]");
		resultMap.put("aliphatic carboxlic acids", "C-C(=O)[O;H1,-]");
		resultMap.put("ketons with exclusion", "[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]");
		resultMap.put("aryl methyl sites", "[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]");
		resultMap.put("carbonyl", "[CX3]=[OX1]");
		resultMap.put("carbonyl with exclusion", "[C!$(C-[OH])]=O	[C!$(C-[OH])]=O");
		resultMap.put("halogens", "[#9,#17,#35,#53]");
		resultMap.put("imidazole rings", "n1cncc1");
		resultMap.put("ether oxygens", "[OD2]([#6])[#6]");
		resultMap.put("tertiary amines", "[NH0,nH0]");
		resultMap.put("methoxy groups", "[OX2](-[#6])-[CH3]");
		resultMap.put("amides", "C(=O)-N");
		resultMap.put("piperzine rings", "N1CCNCC1");
		resultMap.put("aromatic nitrogen", "n");
		resultMap.put("sulfonamides", "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]");
		resultMap.put("benzene rings", "c1ccccc1");
		resultMap.put("aliphatic hydroxyl", "[C!$(C=O)]-[OH]");
		resultMap.put("secondary amines", "[NH1,nH1]");
		resultMap.put("branched alkanes", "[R0;D2][R0;D2][R0;D2][R0;D2]");
		resultMap.put("thiophene rings", "s1cccc1");
		resultMap.put("guanidine groups", "C(=N)(N)N");
		resultMap.put("aliphatic hydroxyl groups with exclusion", "[$(C-[OX2H]);!$([CX3](-[OX2H])=[OX1]);!$([CD4]-[OX2H])]");
		resultMap.put("H-pyrrole nitrogen", "[nH]");
		resultMap.put("aromatic amines", "[nH]");
		resultMap.put("para-hydroxylation sites", "[$([cH]1[cH]cc(c[cH]1)~[$([#8,$([#8]~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)~[$([#7X3,$([#7](~[H,c,C])~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)-!:[$([NX3H,$(NC(=O)[H,c,C])])])]");
		resultMap.put("phenolic OH with exclusion", "[$(c1(-[OX2H])ccccc1);!$(cc-!:[CH2]-[OX2H]);!$(cc-!:C(=O)[O;H1,-]);!$(cc-!:C(=O)-[NH2])]");
		resultMap.put("phenols", "[OX2H]-c1ccccc1");
		resultMap.put("aromatic hydroxyl groups", "c[OH1]");
		resultMap.put("tert-alicyclic amines", "[$([N&R1]1(-C)CCC1),$([N&R1]1(-C)CCCC1),$([N&R1]1(-C)CCCCC1),$([N&R1]1(-C)CCCCCC1),$([N&R1]1(-C)CCCCCCC1)]");
		resultMap.put("XCCNR groups", "[$(N(-[CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)]),$(N(-[CH2][CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)])]");
		resultMap.put("quarternary nitrogens", "[$([NX4+]),$([NX4]=*)]");
		resultMap.put("anilines", "c-[NX3;!$(N=*)]");
		resultMap.put("piperdine rings", "N1CCCCC1");
		resultMap.put("pyridine rings", "n1ccccc1");
		resultMap.put("morpholine rings", "O1CCNCC1");
		resultMap.put("dihydropyridines", "[$([NX3H1]1-C=C-C-C=C1),$([Nv3]1=C-C-C=C-C1),$([Nv3]1=C-C=C-C-C1),$([NX3H1]1-C-C=C-C=C1)]");
		resultMap.put("amidine groups", "C(=N)(-N)-[!#7]");
		resultMap.put("primary amines", "[NH2,nH2]");
		resultMap.put("sulfone groups", "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]");
		resultMap.put("imide groups", "N(-C(=O))-C=O");
		resultMap.put("primary amides", "C(=O)-[NH2]");
		resultMap.put("N functional groups attched to aromatics", "[$(a-[NX3H2]),$(a-[NH1][NH2]),$(a-C(=[OX1])[NH1][NH2]),$(a-C(=[NH])[NH2])]");
		resultMap.put("hydrazone groups", "C=N-[NX3]");
		return resultMap;
	}
}
