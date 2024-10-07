package biotransformer.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.formats.INChIFormat;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.transformation.Biotransformation;
import nu.xom.jaxen.expr.DefaultAbsoluteLocationPath;

public class HandlePolymers {
	public SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	public SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
	HandlePolymerUtils hp_utils = new HandlePolymerUtils();
	HashMap<String, String> aminoAcid_peptide_map;
	//final String smirks_glucose_ish = "[#6]-[#8]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1[#8,#7;A]";
	//final String smirks_glucose_ish = "[#8;A;X2;$([OX2](-C)-C)][#6]-1-[#8]-[#6](-[#6]-*)-[#6](-[#8])-[#6](-[#8])-[#6]-1[#8,#7;A]";
	final String smirks_glucose_ish = "[#8;A;X2;$([OX2]~C)]~C1([H])[#8]C([H])([#6]-[#8])C([H])([#8])C([H])([#8])C1([H])[#8,#7;A]";//"[H]C1([#8;A;X2;$([OX2](-C)-C)])[#8]C([H])([#6]-[#8])C([H])([#8])C([H])([#8])C1([H])[#8,#7;A]";
	final String smirks_gluthion = "[#8]-[#6]-1=[#6]-[#6](-[#6]=O)=[#6]-[#6](-[#8])=[#6]-1-[#8]";
	final String smirks_fructose = "[#8]-[#6]-[#6]-1-[#8]C([#8])([#6]-[#8])[#6](-[#8])-[#6]-1-[#8]";
	//final String smirks_flavonoid = "[#6]-1-[#6]-c2ccccc2-[#8]-[#6]-1-c1ccccc1";
	final String smirks_flavonoid = "[#6]~1~[#6]~c2ccccc2~[#8]~[#6]~1-c1ccccc1";//"[#6]-,=1-[#6]-c2ccccc2-[#8]-[#6]-,=1-c1ccccc1";
	final String smirks_glucose_end = "[H][#8]C1([H])[#8]-[#6](-[#6]-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]1[#8,#7;A]";
	final String smirks_polyphenol = "[#8]-[#6](=O)-c1c(-[#8,#1])c([#8,#1;A])c([#8,#1;A])c([#8,#1;A])c1[#8,#1;A]"; 
	//poly pipetide
	//RNA/DNA

	public static void main(String[] args) throws Exception{
		HandlePolymers hp = new HandlePolymers();
		//String smiles = "OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@H](O)[C@@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O";//HMDB0000163,D-Maltose
		//String smiles = "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O";//HMDB0000258,sucrose
		//String smiles = "OCC1OC(CO)(OC2(O)OC(CO)C(O)C(O)C2O)C(O)C1O";//HMDB0258549,Sucrose alcohol; Not sure what monomer is
		//String smiles = "OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1O";
		//String smiles = "C1=C(C=C(C(=C1O)O)O)C(=O)OC2=CC(=CC(=C2O)O)C(=O)OCC3C(C(C(C(O3)OC(=O)C4=CC(=C(C(=C4)OC(=O)C5=CC(=C(C(=C5)O)O)O)O)O)OC(=O)C6=CC(=C(C(=C6)OC(=O)C7=CC(=C(C(=C7)O)O)O)O)O)OC(=O)C8=CC(=C(C(=C8)OC(=O)C9=CC(=C(C(=C9)O)O)O)O)O)OC(=O)C1=CC(=C(C(=C1)OC(=O)C1=CC(=C(C(=C1)O)O)O)O)O";
		//String smiles = "OCC1OC(CO)(OC2OC(CO)C(OC3OC(CO)C(OC4(CO)OC(CO)C(O)C4O)C(O)C3O)C(O)C2O)C(O)C1O";
		String smiles = "CC(C)C(N)C(=O)NC(C(C)C)C(O)=O";//Peptide dimer
		//String smiles = "N[C@@]([H])(C)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])(CC(C)C)C(=O)O";//peptide
		
		//String smiles = "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)OC5OC(C(O)C(O)C5O)C(O)=O";//"C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O";
		//String smiles = "OC[C@H]1O[C@@H](OC2=COc3cc(O)cc(O)c3C2=O)[C@H](O)[C@@H](O)[C@@H]1O";
		IAtomContainer oneMole = hp.hp_utils.sp.parseSmiles(smiles);
		//AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		System.out.println(hp.isPolymer(oneMole));
		//hp.convertPolysaccaridesToMonoer(oneMole);
		//hp.convertGlucoseMonomers(oneMole);
		ArrayList<String> smirks_list = new ArrayList<>();
		smirks_list.add(hp.smirks_glucose_ish);
		smirks_list.add(hp.smirks_fructose);
		smirks_list.add(hp.smirks_flavonoid);
		smirks_list.add(hp.smirks_glucose_end);
		smirks_list.add(hp.smirks_polyphenol);
		//IAtomContainerSet monomers = hp.hp_utils.convertPolyflavonoidToMonomer(oneMole, hp.smirks_fructose);
		//IAtomContainerSet monomers = hp.hp_utils.convertPolymerToMonomer(oneMole, smirks_list);
		IAtomContainerSet monomers = hp.hp_utils.convertPolymerToMonomer_peptide(oneMole,hp.aminoAcid_peptide_map);
		//IAtomContainerSet monomers = hp.convertPolymerToMonomer(oneMole);
		for(int i = 0; i < monomers.getAtomContainerCount(); i++){
			System.out.println(hp.hp_utils.sg.create(monomers.getAtomContainer(i)));
		}
	}
	
	public HandlePolymers(){
		this.aminoAcid_peptide_map = getAminoAcidsForPeptide();
	}
	/**
	 * Here we only consider if the input molecule contains multiple sugar/flavonoid groups
	 * If the number of sugar (glucose/fructose) gorups + numer of flavonoid groups >=1, we consider it as a polymer and 
	 * will remove the sugar group/flavonoid group till it contains "monomers" only
	 * @param oneMole
	 * @return
	 * @throws Exceptin
	 */
	public boolean isPolymer(IAtomContainer oneMole) throws Exception{	
		if(isPeptide(oneMole)) return true;
		else if(isPolysaccharides(oneMole)) return true;
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		int countGlucose = this.hp_utils.countMonomer(oneMole, this.smirks_glucose_ish);
		int countFructose = this.hp_utils.countMonomer(oneMole, this.smirks_fructose);
		int countFlavonoid = this.hp_utils.countMonomer(oneMole, this.smirks_flavonoid);
		int countpolyphenol = this.hp_utils.countMonomer(oneMole, this.smirks_polyphenol);
		AtomContainerManipulator.suppressHydrogens(oneMole);
		if(this.isPeptide(oneMole)) return true;
		else if(countGlucose + countFructose + countFlavonoid + countpolyphenol>= 2) return true;
		else return false;
	}
	
	/**
	 * Here we only consider if the input molecule contains multiple sugar/flavonoid groups
	 * If the number of sugar (glucose/fructose) gorups + numer of flavonoid groups >=1, we consider it as a polymer and 
	 * will remove the sugar group/flavonoid group till it contains "monomers" only
	 * @param oneMole
	 * @return
	 * @throws Exceptin
	 */
	public boolean isPolymer_moreThanDimer(IAtomContainer oneMole) throws Exception{	
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		int countGlucose = this.hp_utils.countMonomer(oneMole, this.smirks_glucose_ish);
		int countFructose = this.hp_utils.countMonomer(oneMole, this.smirks_fructose);
		int countFlavonoid = this.hp_utils.countMonomer(oneMole, this.smirks_flavonoid);
		AtomContainerManipulator.suppressHydrogens(oneMole);
		if(this.isPeptide(oneMole)) return true;
		else if(countGlucose + countFructose + countFlavonoid >= 3) return true;
		else return false;
	}
	public ArrayList<Biotransformation> converPolymerToMonomer_biotransformation(IAtomContainer molecule) throws Exception{
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = inchiGenFactory.getInChIGenerator(molecule);
		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);	
		molecule.setProperty("InChI", gen.getInchi());
		molecule.setProperty("InChIKey", gen.getInchiKey());
		molecule.setProperty("formula", MolecularFormulaManipulator.getString(formula));
		molecule.setProperty("SMILES", sg.create(molecule));
		molecule.setProperty("Major Isotope Mass", String.valueOf(MolecularFormulaManipulator.getMajorIsotopeMass(formula)));
		//structure.setProperty("RetrievedFromDB", true);
		molecule.setProperty("Reaction ID", String.valueOf("N/A"));
		ALOGPDescriptor aLogpDescriptor = new ALOGPDescriptor(); 
		IDescriptorResult alogp = aLogpDescriptor.calculate(molecule).getValue();
		molecule.setProperty("ALogP", alogp.toString().split(",")[0]);
		IAtomContainer cp = molecule.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(cp);
		IAtomContainerSet substrates = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		substrates.addAtomContainer(molecule);
		IAtomContainerSet monomers = convertPolymerToMonomer(cp);
		
		for(int i = 0; i < monomers.getAtomContainerCount(); i++){
			IAtomContainer structure = monomers.getAtomContainer(i);
			gen = inchiGenFactory.getInChIGenerator(structure);
			formula = MolecularFormulaManipulator.getMolecularFormula(molecule);	
			structure.setProperty("InChI", gen.getInchi());
			structure.setProperty("InChIKey", gen.getInchiKey());
			structure.setProperty("formula", MolecularFormulaManipulator.getString(formula));
			structure.setProperty("SMILES", sg.create(structure));
			structure.setProperty("Major Isotope Mass", String.valueOf(MolecularFormulaManipulator.getMajorIsotopeMass(formula)));
			//structure.setProperty("RetrievedFromDB", true);
			structure.setProperty("Reaction ID", String.valueOf("N/A"));
			aLogpDescriptor = new ALOGPDescriptor(); 
			alogp = aLogpDescriptor.calculate(structure).getValue();
			structure.setProperty("ALogP", alogp.toString().split(",")[0]);
			
		}
		ArrayList<Biotransformation> biotransformations = convertDepolymerizationToBioTransformation(substrates, monomers);
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	/**
	 * This function will take out the monomers with connected atoms from the polymer structure one by one
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet convertPolymerToMonomer(IAtomContainer oneMole) throws Exception{
		//AtomContainerManipulator.suppressHydrogens(oneMole);
		ArrayList<String> backbone_smirksList = new ArrayList<>();
		backbone_smirksList.add(this.smirks_glucose_ish);
		backbone_smirksList.add(this.smirks_fructose);
		backbone_smirksList.add(this.smirks_flavonoid);
		backbone_smirksList.add(this.smirks_glucose_end);
		backbone_smirksList.add(this.smirks_polyphenol);
		IAtomContainerSet result = this.hp_utils.convertPolymerToMonomer(oneMole, backbone_smirksList);
		if(result == null || result.isEmpty()) {
			if(isPolysaccharides(oneMole)) return convertGlucoseMonomers(oneMole);
			else if(isPeptide(oneMole)) return this.hp_utils.convertPolymerToMonomer_peptide(oneMole,aminoAcid_peptide_map);
		}
		return result;
	}
	/**
	 * If the input substrate contains substructures of glucose (also some amino sugar formation) and fructose, and the summation of those functions groups is >=1, then it's considered as polymore
	 * Note that the monomer will not match the two smirks string listed below, because it's -OH group is not -O-C group.
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public boolean isPolysaccharides(IAtomContainer oneMole) throws Exception{
		int totalMatches = 0;
		String smirks_glucose_ish = "[#6]-[#8]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1[#8,#7;A]";
		String smirks_fructose = "[#6]-[#8]C1([#6]-[#8])[#8]-[#6](-[#6]-[#8])-[#6](-[#8])-[#6]1-[#8]";
		int countGlucose = this.hp_utils.countMonomer(oneMole,smirks_glucose_ish);
		int countFructose = this.hp_utils.countMonomer(oneMole, smirks_fructose);
		totalMatches = countGlucose + countFructose;
		if(totalMatches > 1) return true;
		else return false;
	}

	
	public IAtomContainerSet convertGlucoseMonomers(IAtomContainer oneMole) throws Exception{
		String smirks_remove_glucose = "[#6:13]-[#8:3]-[#6:2]-1-[#8:1]-[#6:10](-[#6:11]-[#8:12])-[#6:8](-[#8:9])-[#6:6](-[#8:7])-[#6:4]-1[#8,#7;A:5]>>[H][#8:3]-[#6:2]-1-[#8:1]-[#6:10](-[#6:11]-[#8:12])-[#6:8](-[#8:9])-[#6:6](-[#8:7])-[#6:4]-1[#8,#7;A:5].[H][#8]-[#6:13]";
		IAtomContainerSet finalResult = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet result = getMonomer(oneMole, this.smirks_glucose_ish, smirks_remove_glucose);
		for(int i = 0; i < result.getAtomContainerCount(); i++){
			IAtomContainer recovered = this.sp.parseSmiles(this.sg.create(result.getAtomContainer(i)));
			AtomContainerManipulator.suppressHydrogens(recovered);
			finalResult.addAtomContainer(recovered);
		}
		
		return finalResult;
	}
	/**
	 * This function will keep removing the targetGroup from the query molecule till there is no more matched substructures
	 * @param oneMole
	 * @param smirks_TargetGroup
	 * @param smirks_reaction
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet getMonomer(IAtomContainer oneMole, String smirks_TargetGroup, String smirks_reaction) throws Exception{
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		//If the query molecule is a monomer, then add it into the resultSet and return.
		if(this.hp_utils.countMonomer(oneMole, smirks_TargetGroup) == 0){
			resultSet.addAtomContainer(oneMole);
			return resultSet;		
		}
		else{
			//If there is at least one 
			IAtomContainerSet result = this.hp_utils.removeSubstructure(oneMole, smirks_reaction);
			for(int i = 0; i < result.getAtomContainerCount(); i++){
				IAtomContainer oneMetabolite_temp = result.getAtomContainer(i);
				resultSet.add(getMonomer(oneMetabolite_temp, smirks_reaction, smirks_reaction));
			}
		}
		return resultSet;
	}
	/**
	 * This function will remove the glucose group and fructose group from the query molecule one by one
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet convertPolysaccaridesToMonoer(IAtomContainer oneMole) throws Exception{
		String smirks_remove_glucose = "[#6:13]-[#8:3]-[#6:2]-1-[#8:1]-[#6:10](-[#6:11]-[#8:12])-[#6:8](-[#8:9])-[#6:6](-[#8:7])-[#6:4]-1[#8,#7;A:5]>>[H][#8:3]-[#6:2]-1-[#8:1]-[#6:10](-[#6:11]-[#8:12])-[#6:8](-[#8:9])-[#6:6](-[#8:7])-[#6:4]-1[#8,#7;A:5].[H][#8]-[#6:13]";
		String smirks_remove_fructose = "[#6:13]-[#8:12][C:1]1([#6:10]-[#8:11])[#8:2]-[#6:3](-[#6:8]-[#8:9])-[#6:4](-[#8:7])-[#6:5]1-[#8:6]>>[H][#8:12][C:1]1([#6:10]-[#8:11])[#8:2]-[#6:3](-[#6:8]-[#8:9])-[#6:4](-[#8:7])-[#6:5]1-[#8:6].[H][#8]-[#6:13]";
		IAtomContainerSet finalResult = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet result = this.hp_utils.removeSubstructure(oneMole, smirks_remove_glucose);
		for(int i = 0; i < result.getAtomContainerCount(); i++){
			IAtomContainer oneMetabolite = result.getAtomContainer(i);
			IAtomContainerSet result_2 = this.hp_utils.removeSubstructure(oneMetabolite, smirks_remove_fructose);
			if(result_2!=null && !result_2.isEmpty()){
				for(int j = 0; j < result_2.getAtomContainerCount(); j++){
					finalResult.addAtomContainer(result_2.getAtomContainer(j));
				}
			}
			else finalResult.addAtomContainer(oneMetabolite);
			
		}
		
		return finalResult;
	}
	
	/*
	 * Convert the cyProduct predicted Results to BioTransformations
	 * @param containers
	 * @param nrOfSteps
	 * @param scoreThreshold
	 * @param outputFolder
	 * @param annotate
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> convertDepolymerizationToBioTransformation(IAtomContainerSet substrates, IAtomContainerSet monomers) throws Exception{
		//Biotransformation(IAtomContainerSet substrates, String reactionType, ArrayList<String> enzymeNames, 
		//IAtomContainerSet products, BioSystemName bsysName)
		ArrayList<Biotransformation> convertedBioTrans = new ArrayList<>();
		for(int i = 0; i < monomers.getAtomContainerCount(); i++){
			IAtomContainerSet notIncludedMolecuels = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			IAtomContainer oneMolecule = monomers.getAtomContainer(i);
			notIncludedMolecuels.addAtomContainer(oneMolecule);
			String reactionType = "Depolymerization";
			String enzymeNameString = "Not specified";
			ArrayList<String> enzymeNames = new ArrayList<>();
			enzymeNames.add(enzymeNameString);
			BioSystemName bioSys = BioSystemName.HUMAN;		
			Biotransformation bioTrans = new Biotransformation(substrates, reactionType, enzymeNames, notIncludedMolecuels, bioSys);
			convertedBioTrans.add(bioTrans);
		}
		return convertedBioTrans;
	}
	/**
	 * This function checks if the query molecule is a peptide by:
	 * 1. It must consists of amino acids only
	 * 2. For each of the 20 amino acids, we removed the -OH group is removed from the NCCOOH section and use it SMART string patterns.
	 * All except at most one Hydroxyl group(O) of the molecule should match those patterns.
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public boolean isPeptide(IAtomContainer molecule) throws Exception{
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		ArrayList<Integer> matchedAtomsIdxArray = new ArrayList<>();
		int countAminoAcids = 0;
		for(String name : this.aminoAcid_peptide_map.keySet()){
			String smarts = this.aminoAcid_peptide_map.get(name);
			Pattern pattern = SmartsPattern.create(smarts);
			//int[] matchedAtomIdx = pattern.match(molecule);
			Mappings matchedAtomsMaps = pattern.matchAll(molecule);
			Iterator<int[]> matchedAtomsIterator = matchedAtomsMaps.iterator();
			while(matchedAtomsIterator.hasNext()){
				int[] matchedAtomIdx = matchedAtomsIterator.next();
				boolean newAminoAcid = false;				
				for(int i = 0; i < matchedAtomIdx.length; i++){
					if(!matchedAtomsIdxArray.contains(matchedAtomIdx[i])){
						matchedAtomsIdxArray.add(matchedAtomIdx[i]);
						newAminoAcid = true;					
					}
				}
				if(newAminoAcid){
					countAminoAcids++;
				}
			}	
		}
		int countHydroxyl = 0;
		for(int i = 0; i < molecule.getAtomCount(); i++){
			IAtom oneAtom =  molecule.getAtom(i);
			int index = molecule.indexOf(oneAtom);
			if(!matchedAtomsIdxArray.contains(index)){
				if(oneAtom.getSymbol().equalsIgnoreCase("O")){
					countHydroxyl++;
				}
				if(!oneAtom.getSymbol().equalsIgnoreCase("O") && !oneAtom.getSymbol().equalsIgnoreCase("H"))
					return false;
			}
		}
		if(countHydroxyl>=2 || countAminoAcids <2) return false;
		return true;
	}
	/**
	 * This function creates the HashMap<String,String> for the amino acids patterns used in peptides.
	 * The SMARTS of those amino acids will be used to create monomers from the polymer of peptides
	 * @return
	 */
	public HashMap<String, String> getAminoAcidsForPeptide(){
		HashMap<String, String> resultMap = new HashMap<>();
		resultMap.put("Alanine", "[H]C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Arginine", "[H]\\[#7]=[#6](/[#7]([H])[H])-[#7]([H])C([H])([H])C([H])([H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Asparagine", "[H][#7]([H])-[#6](=O)C([H])([H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Aspartate", "[H][#8]-[#6](=O)C([H])([H])C([H])([#7])[#6]=O");	
		resultMap.put("Cysteine", "[H][#16]C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Glutamine", "[H][#7]([H])-[#6](=O)C([H])([H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Glutamate", "[H][#8]-[#6](=O)C([H])([H])C([H])([H])C([H])([#7])[#6]=O");
		resultMap.put("Glycine", "[H]C([H])([#7])[#6]=O");		
		resultMap.put("Histidine", "[H][#7]-1-[#6]([H])=[#7]-[#6](=[#6]-1[H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Isoleucine", "[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])([#7])[#6]=O");		
		resultMap.put("Leucine", "[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])([#7])[#6]=O");				
		resultMap.put("Lysine", "[H][#7]([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Methionine", "[H]C([H])([H])[#16]C([H])([H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Phenylalanine", "[H][#6]-1=[#6]([H])-[#6]([H])=[#6](-[#6]([H])=[#6]-1[H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Proline", "[H]C1([H])[#7]C([H])([#6]=O)C([H])([H])C1([H])[H]");		
		resultMap.put("Serine", "[H][#8]C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Theonine", "[H][#8]C([H])(C([H])([H])[H])C([H])([#7])[#6]=O");				
		resultMap.put("Tryptophan", "[H][#7]-1-[#6]([H])=[#6](-[#6]-2=[#6]-1-[#6]([H])=[#6]([H])-[#6]([H])=[#6]-2[H])C([H])([H])C([H])([#7])[#6]=O");		
		resultMap.put("Tysosine", "[H][#8]-[#6]-1=[#6]([H])-[#6]([H])=[#6](-[#6]([H])=[#6]-1[H])C([H])([H])C([H])([#7])[#6]=O");				
		resultMap.put("Valine", "[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([#7])[#6]=O");	
		return resultMap;
	}
	
}
