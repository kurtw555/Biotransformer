package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.transformation.Biotransformation;

public class HumanMetabolismHelpingFunctions {
	public SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
	public SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	
	/**
	 * This function will check if the query molecule is a valid gut microbial metabolite.
	 * 
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public boolean isValidGutMetabolite(IAtomContainer oneMole) throws Exception{
		return true;
	}
	
	/**
	 * This function will extract the valid gut microbial biotransformations from all queried ones
	 * @param biotransformations
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> extractValidGutBiotransformations(ArrayList<Biotransformation> biotransformations) throws Exception{
		ArrayList<Biotransformation> result = new ArrayList<>();
		for(int i = 0; i < biotransformations.size(); i++){
			boolean valid_transformation = true;
			Biotransformation biotransformation_one = biotransformations.get(i);
			IAtomContainerSet products_one = extractProductsFromBioTransformation_one(biotransformation_one);
			for(int j = 0; j < products_one.getAtomContainerCount(); j++){
				if(!isValidGutMetabolite(products_one.getAtomContainer(j))){
					valid_transformation = false;
					break;
				}
			}
		
			if(valid_transformation && !result.contains(biotransformation_one)){
				result.add(biotransformation_one);
			}
		}
		
		return result;
	}
	
	/**
	 * This function will extract the unique products for the queried list of biotransformations
	 * Note that if one product occurs in more than one biotransformations, it is stored only once.
	 * @param biotransformations
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet extractProductsFromBiotransformations(ArrayList<Biotransformation> biotransformations) throws Exception{
		IAtomContainerSet result = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<String> storedInChiKey = new ArrayList<>();
		for(Biotransformation oneBiotransformation : biotransformations){
			for(IAtomContainer oneProduct : oneBiotransformation.getProducts().atomContainers()){			
				String inchiKey = oneProduct.getProperty("InChIKey");
				if(!storedInChiKey.contains(inchiKey)){
					storedInChiKey.add(inchiKey);
					result.addAtomContainer(ChemStructureManipulator.preprocessContainer(oneProduct));
				}
			}
		}
		return result;
	}
	
	/**
	 * Extract the unique products from the query biotransformation
	 * @param oneBiotransformation
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet extractProductsFromBioTransformation_one(Biotransformation oneBiotransformation) throws Exception{
		IAtomContainerSet result = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<String> storedInCHiKeyList = new ArrayList<>();
		IAtomContainerSet allProducts = oneBiotransformation.getProducts();
		for(int i = 0; i < allProducts.getAtomContainerCount(); i++){
			IAtomContainer oneProduct = allProducts.getAtomContainer(i);
			String inChiKey = oneProduct.getProperty("InChIKey");
			if(!storedInCHiKeyList.contains(inChiKey)){
				storedInCHiKeyList.add(inChiKey);
				result.addAtomContainer(ChemStructureManipulator.preprocessContainer(oneProduct));
			}
		}
		return result;
	
	}
	
	
	
	
	
	
}
