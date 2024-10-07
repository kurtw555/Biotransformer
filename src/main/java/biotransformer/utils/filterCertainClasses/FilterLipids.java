package biotransformer.utils.filterCertainClasses;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.smsd.ring.RingFinder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.utils.ChemicalClassFinder;

public class FilterLipids {
	
	public void testFA() throws Exception{
		String aliphatic_FA = "O=C(O)CCCCCCC\\C=C/CCCCCC";
		String ring = "[#6]-[#6]-[#6]-1-[#6]-[#6]-[#6]-[#6](-[#6]\\[#6]=[#6]/[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#8])=O)-[#6]-1";
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer molecule = sp.parseSmiles(ring);
		System.out.println(isFattyAcid(molecule));
	}
	
	public boolean isLipid(IAtomContainer molecule) throws Exception{
		if(ChemicalClassFinder.isEtherLipid(molecule)){
			System.out.println("Query molecule is Etherlipid");
		}
		else if(ChemicalClassFinder.isGlyceroLipid(molecule)){
			System.out.println("Query molecule is GlyceroLipid");
		}
		else if(ChemicalClassFinder.isGlycerophosphoLipid(molecule)){
			System.out.println("Query molecule is GlycerophosphoLipid");
		}
		else if(ChemicalClassFinder.isSphingoLipid(molecule)){
			System.out.println("Query molecule is SphingoLipid");
		}
		else if(isFattyAcid(molecule)){
			System.out.println("Query molecule is fatty acid");
		}
		else return false;
		
		return true;
	}
	/**
	 * This function checks if the molecule is a Glycerolipid using two SMIRKS strings:
	 * 1. backbone_smirks_side has the FA chain connects to one side
	 * 2. backbone_smirks_middle has the FA chain connects to the  middle
	 * 
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
//	public boolean isGlycerolipid(IAtomContainer molecule) throws Exception{
//		IAtomContainer cp = molecule.clone();
//		String backbone_smirks_side = "[H]C([#8])([#6]-[#8])[#6]-[#8]-[#6](-[#6])=O";
//		String backbone_smirks_middle = "[H]C([#6]-[#8])([#6]-[#8])[#8]-[#6](-[#6])=O";
//		String backbone_smirks_toRemove = "[H]C([H])([#8])C([H])([#8])C([H])([H])[#8]";
//		Pattern backbone_pattern_side = SmartsPattern.create(backbone_smirks_side);
//		Pattern backbone_pattern_middle = SmartsPattern.create(backbone_smirks_middle);
//		Pattern backbone_pattern_toRemove = SmartsPattern.create(backbone_smirks_toRemove);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(cp);
//		//If this structure doesn't have any side chain connected, then it's not a lipid
//		if(!backbone_pattern_side.matches(cp) && !backbone_pattern_middle.matches(cp)) return false;
//		//If it has the backbone and has at least one FA side chain, then we will remove the backbone and check if the things left contain at least one FA
//		int[] backboneAtomList = backbone_pattern_toRemove.match(cp);
//		IAtom[] atomToRemoveList = new IAtom[backboneAtomList.length];
//		for(int i = 0; i < backboneAtomList.length; i++){
//			int atomIdx = backboneAtomList[i];
//			atomToRemoveList[i] = cp.getAtom(atomIdx);
//		}
//		for(int i = 0; i < atomToRemoveList.length; i++){
//			IAtom atomToRemove = atomToRemoveList[i];
//			if(atomToRemove.getSymbol().equalsIgnoreCase("O")){
//				List<IAtom> neighbor = cp.getConnectedAtomsList(atomToRemove);
//				if(IAtom)
//			}
//			cp.removeAtom(atomToRemoveList[i]);
//		}
//	}
	
	/**
	 * This function defines fatty acids
	 * A fatty acids is:
	 * 1. Contains COOH group at one end
	 * 2. Is aliphatic
	 * 3. should not contains atoms like S, N, P etc.
	 */
	public boolean isFattyAcid(IAtomContainer molecule) throws Exception{		
		AllRingsFinder arf = new AllRingsFinder();
	    IRingSet rs = arf.findAllRings(molecule);
	    //If the structure contains any ring, then it's not fatty acid by definition
	    if(!rs.isEmpty()) return false;
	    for(int i = 0; i < molecule.getAtomCount(); i++){
	    	IAtom oneAtom = molecule.getAtom(i);
	    	//If it contains any 
	    	if(!oneAtom.getSymbol().equalsIgnoreCase("H") && oneAtom.getSymbol().equalsIgnoreCase("C") && oneAtom.getSymbol().equalsIgnoreCase("O")){
	    		return false;
	    	}
	    }
	    //If the structure doesn't contain any ring, then we check if there is a COOH
	    String smirks_COOH = "[H][#8]-[#6]=O";
	    String smirks_OCO = "[#6]-[#8]-[#6]";
	    String smirks_OO = "[#8]-[#8]";
	    Pattern pattern_COOH = SmartsPattern.create(smirks_COOH);
	    Pattern pattern_OCO = SmartsPattern.create(smirks_OCO);
	    Pattern pattern_OO = SmartsPattern.create(smirks_OO);
	    AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
	    //If the molecule doesn't have COOH, then it's not fatty acid
	    if(!pattern_COOH.matches(molecule)) return false;
	    //If the molecule contains ester or OO, then return again it's not a faaty acid
	    if(pattern_OCO.matches(molecule) || pattern_OO.matches(molecule)) return false;
	    //If the molecule passed all tests above, then it's a fatty acid
	    AtomContainerManipulator.suppressHydrogens(molecule);
	    return true;
	}
	
	
	
	
	
	
}
