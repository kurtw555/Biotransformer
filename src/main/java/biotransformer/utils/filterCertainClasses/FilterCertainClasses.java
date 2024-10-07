package biotransformer.utils.filterCertainClasses;

import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.dbrelevant.RetriveFromDB;
import reactantpredictor.ReactantPred;

public class FilterCertainClasses {
	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	FilterAminoAcids filterAminoAcids = new FilterAminoAcids();
	FilterCarbohydrates filterCarbohydrates = new FilterCarbohydrates();
	FilterLipids filterLipids = new FilterLipids();	
	
	/**
	 * This function checks if the query molecule should go to retrieve metabolites from Pathbank directly
	 * If yes, then it will goes to retrieving from pathbank process -- itself and its metabolites will not be predicted by Biotransformer again
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public boolean retriveFromDB(IAtomContainer molecule) throws Exception{
		InChIGeneratorFactory inChIGenerator = InChIGeneratorFactory.getInstance();
		InChIGenerator inchiGen = inChIGenerator.getInChIGenerator(molecule);
		String inChIKey = inchiGen.getInchiKey();
		RetriveFromDB rfpb = new RetriveFromDB("hmdb", "all", false, 30, false);
		//We add a manual fix for D-Glucose
		if(inChIKey.equals("WQZGKKKJIJFFOK-GASJEMHNSA-N") || inChIKey.equals("WQZGKKKJIJFFOK-DVKNGEFBSA-N") || inChIKey.equals("WQZGKKKJIJFFOK-UHFFFAOYNA-N") || inChIKey.equals("WQZGKKKJIJFFOK-UHFFFAOYSA-N")) inChIKey = "WQZGKKKJIJFFOK-VFUOTHLCSA-N";
		return rfpb.isEndogenous(inChIKey);
//		return true;
//		if(this.filterAminoAcids.isAminoacid(molecule)) return true;		
//		else if(this.filterCarbohydrates.isCarbohydrates(molecule)) return true;
//		else if(this.filterLipids.isLipid(molecule)) return true;
//		else return false;
	}
	
	public boolean retriveFromPathbank(String smiles) throws Exception{
		IAtomContainer molecule = this.sp.parseSmiles(smiles);
		return retriveFromDB(molecule);
	}
	/**
	 * This function will partition the query set of molecules.
	 * ArrayList<> = {SetOfMoleculesSatisfyDBCondition, SetOfMoleculesNoStatisfyDBCondition}
	 * if there is no molecule in resultSet_db, then resultArray.get(0) will be empty and that notifies the other funcitons that there are no matched candidates from DB
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public ArrayList<IAtomContainerSet> getMoleculesSatifiesDBCondition(IAtomContainerSet molecules) throws Exception{
		ArrayList<IAtomContainerSet> resultArray = new ArrayList<>();
		IAtomContainerSet resultSet_db = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet resultSet_no_db = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer molecule = molecules.getAtomContainer(i);
			if(retriveFromDB(molecule)){
				resultSet_db.addAtomContainer(molecule);
			}
			else resultSet_no_db.addAtomContainer(molecule);
		}
		resultArray.add(resultSet_db);
		resultArray.add(resultSet_no_db);
		return resultArray;
	}
}
