package biotransformer.utils.filterCertainClasses;

import java.util.HashMap;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public class FilterCarbohydrates {
	/**
	 * This function checks if the molecule's formula matches Cm(H2O)n. If not, then definitely it's not a carbohydrate
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public boolean isCarbohydrates(IAtomContainer molecule) throws Exception{
		//IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);	
		//String formula_string = MolecularFormulaManipulator.getString(formula);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		Double countCarbon = 0.0;
		Double countOxygen = 0.0;
		Double countHydrogen = 0.0;
		for(int i = 0; i < molecule.getAtomCount(); i++){
			IAtom oneAtom = molecule.getAtom(i);
			String symbol = oneAtom.getSymbol();
			if(symbol.equalsIgnoreCase("C")) countCarbon++;
			else if(symbol.equalsIgnoreCase("O")) countOxygen++;
			else if(symbol.equalsIgnoreCase("H")) countHydrogen++;
			else return false;
			
		}
		if(countHydrogen/countOxygen == 2.0){
			return true;
		}
		else return false;
	}
	

	/**
	 * This function checks if the molecule is made of monomer sugars. If yes, then it's a  carbohydrate, otherwise not.
	 * This function will only be called after the molecule pass checkFormula_Carbohydrates check which is easier and faster.
	 * Because 1. none of the SMIRKS used is a subset of another one; 2. when the C-O-C connection is formed, a water is removed,
	 * assume there are m matches for glucose, n matches for fructose and k matches for xylose, we will check the number of atoms = 12*m + 12*n + 10*k - (m+n+k - 1).
	 * If that number = numer of non-hydrogen atoms in the moelcule, then return true; otherwise return false;
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
//	public boolean checkMonomer_Carbohydrates(IAtomContainer molecule) throws Exception{
//		//
//		String glucose_ring_SMIRKS = "[#8]-[#6]-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]";
//		String fructose_ring_SMIRKS = "[#8]-[#6]-[#6]-1-[#8]C([#8])([#6]-[#8])[#6](-[#8])-[#6]-1-[#8]";
//		String xylose_ring_SMIRKS = "[#8]-[#6]-1-[#6]-[#8]-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]";
//		return true;
//		
//	}
	
}
