package biotransformer.validateModels;

import java.util.ArrayList;
import java.util.HashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class IsEndProductSMARTS {
	public HashMap<String, String> endProductSubstructuresMap = new HashMap<>();
	
	public IsEndProductSMARTS(){
		this.endProductSubstructuresMap = getEndProductStructureMap();
	}
	
	
	public static void main(String[] args) throws Exception{
		IsEndProductSMARTS isEndProduct = new IsEndProductSMARTS();
		//String smiles = "C1(C(C2=C(C=C(C=C2OC1C3=CC(=C(C=C3)OC4OC(C(O)C(O)C4O)C(O)=O)O)O)OC5OC(C(O)C(O)C5O)C(O)=O)=O)O";
		//String smiles = "OC1C(O)C(OC2=C(O)C=C(C=C2)C2OC3=CC(O)=CC(O)=C3C(=O)C2O)OC(C1O)C(O)=O";
		//String smiles = "C1(C(C2=C(C=C(C=C2OC1C3=CC(=C(C=C3)O)OS(O)(=O)=O)O)OC4OC(C(O)C(O)C4O)C(O)=O)=O)O";
		String smiles = "C(=O)(O)C1=CC(=CC(C1O)SCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O)O";
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer oneMole = sp.parseSmiles(smiles);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		System.out.println(isEndProduct.isEndProduct(oneMole));
	}
	/**
	 * If it's not a phase II metabolism metabolite.
	 * If the molecule contains at least two endProductSubstructres, then we assume it is an end product
	 * Note:
	 * 1. If there is a Glutathione, then the glycine will also be matched. In this case, the counter will be >=2. However, we think Glutathione should be the end of the metabolism, so it's fine.
	 * 
	 * 
	 * @return
	 */
	public boolean isEndProduct(IAtomContainer oneMole) throws Exception{
		if(oneMole.getProperty("isEndProduct") != null) return true;
		int counter = 0;
		//In some cases, these conjugation metabolism can occur multiple times, say 2.
		//Hence we set the total number threshold as 2, an this may work for human metabolism
		for(String name : this.endProductSubstructuresMap.keySet()){
			String smirks = this.endProductSubstructuresMap.get(name);
			Pattern pattern = SmartsPattern.create(smirks);
			Mappings mappings = pattern.matchAll(oneMole);
			if(name.equals("Glucuronidation")||name.equals("Sulfation")||name.equals("Phosphorylation")||name.equals("Glutathione")){
				if(countPossibleUniqueMappings(mappings) > 2){
					return true;
				}
			}
			counter += countPossibleUniqueMappings(mappings);
		}
		if(counter >= 2) return true;
		else return false;
	}
	
	public int countPossibleUniqueMappings(Mappings mappings){
		int[][] matrix = mappings.toArray();
		int counter = 0;
		ArrayList<ArrayList<Integer>> checkArray = new ArrayList<>();
		for(int i = 0; i < matrix.length; i++){
			int[] oneRow = matrix[i];
			boolean duplicate = false;
			ArrayList<Integer> temp = new ArrayList<>();
			if(checkArray.isEmpty()){
				ArrayList<Integer> oneArray = new ArrayList<>();
				for(int j = 0; j < oneRow.length; j++){
					int oneEntry = oneRow[j];
					oneArray.add(oneEntry);
				}
				counter++;
				checkArray.add(oneArray);
			}
			else{
				for(int t = 0; t < checkArray.size(); t++){
					int matchedAtoms = 0;
					ArrayList<Integer> oneArray = checkArray.get(t);					
					for(int j = 0; j < oneRow.length; j++){
						if(t==0) temp.add(oneRow[j]);//we only want to create the temp array for the current row once
						if(oneArray.contains(oneRow[j])) matchedAtoms++;
					}
					//If more than half of the atoms are duplicated, then that row is not valid -- it's the same as the previous one
					if(matchedAtoms >= oneRow.length/2){
						continue;
					}
					else{
						duplicate = true;
						break;
					}					
				}
				if(duplicate) counter++;
				if(!checkArray.contains(temp)) checkArray.add(temp);
			}			
		}
		return counter;
	}
	/**
	 * This function checks if the generated metabolite contains some obviously unreasonable substructure
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public boolean containEndSubstructure(IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		boolean contains = false;
		for(String key : this.endProductSubstructuresMap.keySet()){
			String smarts = this.endProductSubstructuresMap.get(key);
			Pattern pattern = SmartsPattern.create(smarts);
			if(pattern.matches(oneMole)){
				contains = true;
				break;
			}
			
		}
		AtomContainerManipulator.suppressHydrogens(oneMole);
		return contains;
	}
	
	public HashMap<String, String> getEndProductStructureMap(){
		HashMap<String, String> resultMap = new HashMap<>();
		//resultMap.put("Glucuronidation", "[H][#8]-[#6](=O)-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]");
		resultMap.put("Glucuronidation", "[H][#8]-[#6](=O)-[#6]-1-[#8]-[#6](-*)-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8]");
		resultMap.put("Sulfation", "[H][#8]S(=O)(=O)[#8;X2]-*");
		resultMap.put("Phosphorylation", "[#8;A;X2]P([#8])([#8])=O");
		resultMap.put("Carnitine", "[#6]-[#6](=O)-[#8]-[#6](-[#6]-[#6](-[#8])=O)-[#6][N+]([#6])([#6])[#6]");		
		resultMap.put("Glycine", "[H][#8]-[#6](=O)-[#6]-[#7]-[#6](-[#6])=O");
		resultMap.put("Taurine", "[#6]-[#6](=O)-[#7]-[#6]-[#6]S([#8])(=O)=O");
		resultMap.put("Glutathione", "[H][#8]-[#6](=O)-[#6]-[#7]([H])-[#6](=O)-[#6](-[#6]-[#16])-[#7]([H])-[#6](=O)-[#6]-[#6]-[#6](-[#7]([H])[H])-[#6](=O)-[#8][H]");
		return resultMap;
	}
}
