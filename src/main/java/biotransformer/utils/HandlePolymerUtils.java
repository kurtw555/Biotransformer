package biotransformer.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.query.SmartsPatternAmbit;
import biotransformer.utils.filterCertainClasses.FilterAminoAcids;

public class HandlePolymerUtils {
	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
	protected SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
	PathTools pt = new PathTools();
	/**
	 * This function will check the number of unique matches for the query pattern which is described by the input smirks_string
	 * @param oneMole
	 * @param smirks_string
	 * @return
	 * @throws Exception
	 */
	public int countMonomer(IAtomContainer oneMole, String smirks_string) throws Exception{
		Pattern smirks_pattern = SmartsPattern.create(smirks_string);
		Mappings pattern_mappings = smirks_pattern.matchAll(oneMole);
		//int[][] check = pattern_mappings.toArray();
		return pattern_mappings.countUnique();
	}

	/**
	 * This function will remove the query structure specified by the smirks_string from the query molecule
	 * It then returns the result structures.
	 * @param oneMole
	 * @param smirks_string
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet removeSubstructure(IAtomContainer oneMole, String smirks_string) throws Exception{
		ArrayList<String> checkExist = new ArrayList<>();
		IAtomContainerSet resultStructures = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		SMIRKSReaction oneReaction =  this.smrkMan.parse(smirks_string);
		IAtomContainerSet oneMetabolites_oneSet = this.smrkMan.applyTransformationWithSingleCopyForEachPos(oneMole, null, oneReaction);
		if(oneMetabolites_oneSet == null) return null;
		for(int i = 0; i < oneMetabolites_oneSet.getAtomContainerCount(); i++){
			IAtomContainer oneMetabolites_one = oneMetabolites_oneSet.getAtomContainer(i);
			IAtomContainerSet metabolites = ConnectivityChecker.partitionIntoMolecules(oneMetabolites_one);
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				String oneMetabolite_smiles = this.sg.create(metabolites.getAtomContainer(j));
				if(!checkExist.contains(oneMetabolite_smiles)){
					resultStructures.addAtomContainer(metabolites.getAtomContainer(j));
					checkExist.add(oneMetabolite_smiles);
				}
				
			}
		}
		return resultStructures;
	}
	/**
	 * Remember to remove targetStructureIdxList from the allMaps when generate the substructure for the targetStructureIdxList
	 * AllMaps doesn't contain the targetStructureIdxList
	 * We are build the structure using targetStructureIdxList and all atoms that are not in the other mapped backbones
	 * @param molecule
	 * @param targetStructureIdxList
	 * @param allMaps
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet formStrucutreUsingIdx(IAtomContainer molecule, ArrayList<Integer> targetStructureIdxList, ArrayList<ArrayList<Integer>> allMaps) throws Exception{
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<Integer> allAtomIdxInMonomers = new ArrayList<>();
		for(int i = 0; i < allMaps.size(); i++){
			for(int j = 0; j < allMaps.get(i).size(); j++){
				if(!allAtomIdxInMonomers.contains(allMaps.get(i).get(j))) allAtomIdxInMonomers.add(allMaps.get(i).get(j));				
			}
		}
		IAtomContainer oneMole = molecule.clone();		
		IAtomContainer temp = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		//Add the backbone's atoms into the substructure
		for(int i = 0; i < targetStructureIdxList.size(); i++){			
			IAtom theAtom = oneMole.getAtom(targetStructureIdxList.get(i));
			if(!temp.contains(theAtom)) temp.addAtom(theAtom);
		}
		//Add the atoms that are connected to the atom in the substructure and not in any backbones into the temp
		for(int i = 0; i < oneMole.getAtomCount(); i++){			
			IAtom oneAtom = oneMole.getAtom(i);			
			Integer idx = oneMole.indexOf(oneAtom);
			if(oneAtom.getSymbol().equalsIgnoreCase("H") && !temp.contains(oneAtom)) {
				if(temp.contains(oneMole.getConnectedAtomsList(oneAtom).get(0))) temp.addAtom(oneAtom);
			}
			else if(!allAtomIdxInMonomers.contains(idx) && !temp.contains(oneAtom)){
				if(hasGoodPath(oneMole, temp, idx, allAtomIdxInMonomers)){
					temp.addAtom(oneAtom);
				}
				else continue;
				//Handle ester group
				if(oneAtom.getSymbol().equalsIgnoreCase("S") || oneAtom.getSymbol().equalsIgnoreCase("C") || oneAtom.getSymbol().equalsIgnoreCase("P")) {
					List<IAtom> checkEster = oneMole.getConnectedAtomsList(oneAtom);
					boolean isKetone = false;
					for(int j = 0; j < checkEster.size(); j++) {
						if(checkEster.get(j).getSymbol().equalsIgnoreCase("O")) {
							IAtom checkOxygen =checkEster.get(j);
							if(oneMole.getBond(oneAtom, checkOxygen).getOrder() == IBond.Order.DOUBLE) {
								isKetone = true;
								break;
							}
						}						
					}
					if(isKetone) {
						for(int j = 0; j < checkEster.size(); j++) {
							if(checkEster.get(j).getSymbol().equalsIgnoreCase("O")) {
								IAtom checkOxygen =checkEster.get(j);
								if(oneMole.getBond(oneAtom, checkOxygen).getOrder() == IBond.Order.SINGLE && allAtomIdxInMonomers.contains(oneMole.indexOf(checkOxygen))) {
									IAtom fakeAtom = new Atom("O");
									//fakeAtom.setImplicitHydrogenCount(1);
									IBond fakeBond = new Bond(oneAtom, fakeAtom, oneMole.getBond(oneAtom, checkOxygen).getOrder());
									temp.addAtom(fakeAtom);
									temp.addBond(fakeBond);
									//temp.addAtom(checkOxygen);
									break;
								}
							}						
						}
					}
				}
			}
		}
		
		for(int i = 0; i < temp.getAtomCount(); i++){
			for(int j = 0; j < temp.getAtomCount(); j++){
				IBond oneBond = oneMole.getBond(temp.getAtom(i), temp.getAtom(j));
				if(oneBond != null && !temp.contains(oneBond)) temp.addBond(oneBond);
			}
		}
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(temp);
		for(int i = 0; i < temp.getAtomCount(); i++){
			IAtom oneAtom = temp.getAtom(i);		
			if(oneAtom.getSymbol().equals("H")) continue;
			List<IAtom> neighAtoms = temp.getConnectedAtomsList(oneAtom);
			Double bondOrderSum = (double) oneAtom.getImplicitHydrogenCount();
			for(int j = 0; j < neighAtoms.size(); j++){
				IBond theBond = temp.getBond(oneAtom, neighAtoms.get(j));
				bondOrderSum = bondOrderSum + theBond.getOrder().numeric();
			}
			
			//For some reason, sometimes this occurs. The valency is null. (e.g. surfur). We want to keep adding Hydrogens, but not infinite hydrogens. Hence we set max H = 2 here
			int counter = 0;
			while(oneAtom.getValency() == null && counter < 3) {
				IAtom fakeAtom = new Atom("H");
				IBond fakeBond = new Bond(oneAtom, fakeAtom, IBond.Order.SINGLE);
				temp.addAtom(fakeAtom);
				temp.addBond(fakeBond);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(temp);
				bondOrderSum = bondOrderSum + fakeBond.getOrder().numeric();
				counter++;
			}
			int valence = oneAtom.getValency();
			if(valence > bondOrderSum + oneAtom.getFormalCharge()){
				int diffHydrogen = (int) (valence - bondOrderSum - oneAtom.getFormalCharge());

				oneAtom.setImplicitHydrogenCount(diffHydrogen);
			}
		}
		ArrayList<Integer> rest = new ArrayList<>();
		for(int i = 0; i < oneMole.getAtomCount(); i++) {
			IAtom oneAtom = oneMole.getAtom(i);
			if(!temp.contains(oneAtom)) {
				rest.add(oneMole.indexOf(oneAtom));
				List<IAtom> neighborList = oneMole.getConnectedAtomsList(oneAtom);
				for(int j = 0; j < neighborList.size(); j++) {
					if(neighborList.get(j).getSymbol().equals("O") && temp.contains(neighborList.get(j))){
						rest.add(oneMole.indexOf(neighborList.get(j)));
					}
				}
			}
		}
		//We don't want mixture. The monomer should be removed from the end
		if(temp != null && !temp.isEmpty() && !this.sg.create(temp).contains(".")) {		
			resultSet.addAtomContainer(temp);
		}		
		IAtomContainer rest_struct = createMonomerForBackbone(rest, oneMole);
		if(rest_struct != null && !rest_struct.isEmpty() && !this.sg.create(rest_struct).contains(".")) {
			resultSet.addAtomContainer(rest_struct);
		}
		return resultSet;
	}
	
	/**
	 * This function will check if there is a path from any atom in the temp to the targetAtom.
	 * If there is a valid one, then the targetAtom can be added into the temp
	 * Note that to be a goodPath, none of the atoms on the path can be in the ArrayList<Integer> allAtomIdxInMonomers
	 * @param oneMole
	 * @param temp
	 * @param targetAtomIdx
	 * @param allAtomIdxInMonomers
	 * @return
	 * @throws Exception
	 */
	public boolean hasGoodPath(IAtomContainer oneMole, IAtomContainer temp, Integer targetAtomIdx, ArrayList<Integer> allAtomIdxInMonomers) throws Exception{
		for(int i = 0; i < temp.getAtomCount(); i++){
			IAtom oneAtom = temp.getAtom(i);
			if(!oneMole.contains(oneAtom)) continue;
			IAtom targetAtom = oneMole.getAtom(targetAtomIdx);
			List<List<IAtom>> allPaths = this.pt.getAllPaths(oneMole, oneAtom, targetAtom);
			boolean isValidPath = false;
			for(int j = 0; j < allPaths.size(); j++){
				boolean isValidPath_temp = true;
				for(int k = 0; k < allPaths.get(j).size(); k++){
					Integer onePathAtomIdx = oneMole.indexOf(allPaths.get(j).get(k));
					if(allAtomIdxInMonomers.contains(onePathAtomIdx)){
						isValidPath_temp = false;
						break;
					}
				}
				//If found one valid path from the tempAtom in the temp to the targetAtom
				if(isValidPath_temp){
					isValidPath = true;
					break;
				}
			}
			if(isValidPath){
				return true;
			}
		}
		return false;				
	}
	
	/**
	 * This function will convert a polymer to monomers using the input ArrayList<String> backbone_smirksList
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet convertPolymerToMonomer(IAtomContainer oneMole, ArrayList<String> backbone_smirksList) throws Exception{
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		ArrayList<String> existedInChiKeyList = new ArrayList<>();
		//String backbone_smirks = "[#6]-1-[#6]-c2ccccc2-[#8]-[#6]-1-c1ccccc1";
		//SMIRKSReaction oneReaction =  this.hp_utils.smrkMan.parse(backbone_smirks);
		ArrayList<ArrayList<Integer>> atomIdxList = new ArrayList<>();
		for(int i = 0; i < backbone_smirksList.size(); i++){
			String backbond_smirks = backbone_smirksList.get(i);
			Pattern smirks_pattern = SmartsPattern.create(backbond_smirks);
			Mappings pattern_mappings = smirks_pattern.matchAll(oneMole);										
			for(int[] oneMap : pattern_mappings.uniqueAtoms()){
				ArrayList<Integer> oneList = new ArrayList<>();
				for(int j = 0; j < oneMap.length; j++){
					oneList.add(oneMap[j]);
				}
				atomIdxList.add(oneList);
			}
		}
		for(int i = 0; i < atomIdxList.size(); i++){
			ArrayList<Integer> oneTargetList = atomIdxList.get(i);
			ArrayList<ArrayList<Integer>> temp = (ArrayList<ArrayList<Integer>>) atomIdxList.clone();
			temp.remove(i);
			IAtomContainerSet monomers = formStrucutreUsingIdx(oneMole, oneTargetList, temp);
			for(int j = 0; j < monomers.getAtomContainerCount(); j++) {			
				IAtomContainer monomer = monomers.getAtomContainer(j);
				AtomContainerManipulator.suppressHydrogens(monomer);
				InChIGenerator gen = inchiGenFactory.getInChIGenerator(monomer);
				String inChiKey = gen.getInchiKey();
				if(!existedInChiKeyList.contains(inChiKey)){
					resultSet.addAtomContainer(monomer);
					existedInChiKeyList.add(inChiKey);
				}
			}
		}
//		testFormSubstructure(oneMole, atomIdxList.get(0));
		return resultSet;
	}
	/**
	 * This function will handles peptides only. 
	 * It will extract monomers from the peptids. Please note that:
	 * 1. Each monomer must be one of the 20 amino acids.
	 * 2. No matter how many amino acids A are in the peptide, only one A will be returned as one monomer.
	 * @param oneMole
	 * @param amino_patternList
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet convertPolymerToMonomer_peptide(IAtomContainer oneMole_query, HashMap<String, String> amino_patternList) throws Exception{
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainer oneMole = oneMole_query.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		for(String name : amino_patternList.keySet()){
			String smarts = amino_patternList.get(name);
			SmartsPattern pattern = SmartsPattern.create(smarts);
			Mappings all_maps = pattern.matchAll(oneMole);
			int[] oneMap;
			if(all_maps != null && all_maps.iterator().hasNext()){
				oneMap = all_maps.iterator().next();
			}
			else continue;
			ArrayList<Integer> listOfAtoms = new ArrayList<>();
			for(int i = 0; i < oneMap.length; i++){
				IAtom oneAtom = oneMole.getAtom(oneMap[i]);
				if(!listOfAtoms.contains(oneAtom)) listOfAtoms.add(oneMole.indexOf(oneAtom));
			}
			
			IAtomContainer monomer = createAminoAcidFromPeptide(listOfAtoms, oneMole.clone());
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(monomer);			
			IAtomContainer recovered =  this.sp.parseSmiles(this.sg.create(monomer));
			AtomContainerManipulator.suppressHydrogens(recovered);
			resultSet.addAtomContainer(recovered);
		}
		return resultSet;
	}
	/**
	 * This function will return the recovered amino acid monomer from the query peptide (oneMole) using the matched atoms stored in the ArrayList<IAtom> atomList
	 * @param atomList
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public IAtomContainer createAminoAcidFromPeptide(ArrayList<Integer> atomList, IAtomContainer oneMole) throws Exception{
		IAtomContainer amino_struct = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		String key_pattern_smirks = "[#7]-[#6]-[#6]=O";
		for(int i = 0; i < atomList.size(); i++){
			if(!amino_struct.contains(oneMole.getAtom(atomList.get(i)))) amino_struct.addAtom(oneMole.getAtom(atomList.get(i)));
		}
		for(int i = 0; i < amino_struct.getAtomCount(); i++){
			IAtom atom_one = amino_struct.getAtom(i);
			for(int j = 0; j < amino_struct.getAtomCount(); j++){
				IAtom atom_two = amino_struct.getAtom(j);
				IBond oneBond = oneMole.getBond(atom_one, atom_two);
				if(oneBond!=null && !amino_struct.contains(oneBond)) amino_struct.addBond(oneBond);
			}
		}
		SmartsPattern pattern = SmartsPattern.create(key_pattern_smirks);
		Mappings all_maps = pattern.matchAll(amino_struct);
		if(all_maps == null) throw new Exception("Cannot find the key patterns in the query peptide fragments");
		Iterator<int[]> iterator = all_maps.iterator();
		while(iterator.hasNext()){
			int[] oneMap = iterator.next();
			for(int i = 0; i < oneMap.length; i++){
				IAtom oneAtom = amino_struct.getAtom(oneMap[i]);
				int countImplicitHydrogen = 0;
				if(oneAtom.getSymbol().equalsIgnoreCase("N")){
					List<IAtom> neighbors = oneMole.getConnectedAtomsList(oneAtom);
					for(int j = 0; j < neighbors.size(); j++){
						IAtom checkAtom = neighbors.get(j);
						//If there is a neighbor atom that is not in the amino structure but in the original peptide molecule,
						//then it means it should be replaced with Hydrogen to recover to the original state from peptide bond
						if(!amino_struct.contains(checkAtom) && oneMole.contains(checkAtom)){
							countImplicitHydrogen++;							
						}
					}
					oneAtom.setImplicitHydrogenCount(countImplicitHydrogen);
				}
				if(oneAtom.getSymbol().equalsIgnoreCase("C")){
					List<IAtom> neighbors = oneMole.getConnectedAtomsList(oneAtom);
					for(int j = 0; j < neighbors.size(); j++){
						IAtom checkAtom = neighbors.get(j);
						//If there is a neighbor atom that is not in the amino structure but in the original peptide molecule,
						//then it means it should be replaced with Hydrogen to recover to the original state from peptide bond
						if(!amino_struct.contains(checkAtom) && oneMole.contains(checkAtom)){
							IBond theBond = oneMole.getBond(oneAtom, checkAtom);
							IAtom fakeAtom = new Atom("O");
							fakeAtom.setImplicitHydrogenCount(1);						
							IBond fakeBond = new Bond(oneAtom, fakeAtom, theBond.getOrder());
							amino_struct.addAtom(fakeAtom);
							amino_struct.addBond(fakeBond);
						}
					}
				}
			}
		}
		return amino_struct;
	}
	
	public IAtomContainer createMonomerForBackbone(ArrayList<Integer> oneTargetList, IAtomContainer oneMole_origin) throws Exception{
		IAtomContainer monomer = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		IAtomContainer oneMole = oneMole_origin;
		//Add all atoms in the backbone into the container
		for(int i = 0; i < oneTargetList.size(); i++) {
			IAtom oneAtom = oneMole.getAtom(oneTargetList.get(i));
			if(!monomer.contains(oneAtom)) monomer.addAtom(oneAtom);
		}
		//Add all bonds in the backbone into the container
		for(int i = 0; i < monomer.getAtomCount(); i++) {
			for(int j = 0; j < monomer.getAtomCount(); j++) {
				if(oneMole.getBond(monomer.getAtom(i), monomer.getAtom(j)) != null && !monomer.contains(oneMole.getBond(monomer.getAtom(i), monomer.getAtom(j)))) {
					monomer.addBond(oneMole.getBond(monomer.getAtom(i), monomer.getAtom(j)));
				}
			}
		}
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(monomer);
		//Clean up the atoms connect to the rest of the original molecule
		for(int i = 0; i < monomer.getAtomCount(); i++){
			IAtom oneAtom = monomer.getAtom(i);		
			if(oneAtom.getSymbol().equals("H")) continue;
			List<IAtom> neighAtoms = monomer.getConnectedAtomsList(oneAtom);
			Double bondOrderSum = (double) oneAtom.getImplicitHydrogenCount();
			for(int j = 0; j < neighAtoms.size(); j++){
				IBond theBond = monomer.getBond(oneAtom, neighAtoms.get(j));
				bondOrderSum = bondOrderSum + theBond.getOrder().numeric();
			}
			//For some reason, sometimes this occurs. The valency is null. (e.g. surfur). We want to keep adding Hydrogens, but not infinite hydrogens. Hence we set max H = 2 here
			if(oneAtom.getValency() == null) {
				oneAtom.setImplicitHydrogenCount(1);
				//AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(monomer);
				AtomValenceDescriptor aValenceDes = new AtomValenceDescriptor();
				DescriptorValue valence_calculated = aValenceDes.calculate(oneAtom, monomer);
				Integer vl = Integer.parseInt(valence_calculated.getValue().toString());
				oneAtom.setValency(vl);
			}						
			int valence = oneAtom.getValency();
			if(valence > bondOrderSum + oneAtom.getFormalCharge()){
				int diffHydrogen = (int) (valence - bondOrderSum - oneAtom.getFormalCharge());
				oneAtom.setImplicitHydrogenCount(diffHydrogen);
			}
		}
		return monomer;
	}
	
	
	
	
}
