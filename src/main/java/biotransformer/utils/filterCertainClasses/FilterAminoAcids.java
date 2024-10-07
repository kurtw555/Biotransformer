package biotransformer.utils.filterCertainClasses;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.bcel.generic.NEWARRAY;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class FilterAminoAcids {
	HashMap<String, String> aminoAcidMap;
	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	public static void main(String[] args) throws Exception{
		FilterAminoAcids faa = new FilterAminoAcids();		
		String smiles = "CC(C)C(N)C(O)=O";
		//String smiles = "CC(C)C(N)C(=O)NC(C(C)C)C(O)=O";//Dimer
		//String smiles = "N[C@@]([H])(C)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])(CC(C)C)C(=O)O";//peptide
		//String smiles = "[H][#6](-[#7])-[#6](=O)-[#7][C@@]([H])([#6]-[#16])[#6](=O)-[#7][C@]([H])([#6](=O)-[#7][C@@]([H])([#6]-[#6](-[#6])-[#6])[#6](-[#8])=O)[C@@]([H])([#6])[#8]";//not peptide, with a none of 20 amino acids
		//String smiles = "[H][C@]([#6])([#8])[C@]([H])([#7]-[#6](=O)[C@]([H])([#6]-[#16])[#7]-[#6](=O)[C@@]([H])([#7])[#6]-[#6])[#6](=O)-[#7][C@@]([H])([#6]-[#6](-[#6])-[#6])[#6](-[#8])=O";//Not peptide.with extra branch
		IAtomContainer molecule = faa.sp.parseSmiles(smiles);
		//System.out.println(faa.isAminoacid(molecule));
	}
	public FilterAminoAcids(){
		this.aminoAcidMap = getAminoAcidMap();
		//this.aminoAcid_peptide_map = getAminoAcidsForPeptide();
	}
	
	public boolean isAminoacid(IAtomContainer molecule) throws Exception{
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		for(String name : this.aminoAcidMap.keySet()){
			String smarts = this.aminoAcidMap.get(name);
			Pattern pattern = SmartsPattern.create(smarts);
			int[] matchedAtomIdx = pattern.match(molecule);
			if(matchedAtomIdx.length == molecule.getAtomCount()){
				System.out.println("Matched " + name);
				return true;
			}
			
		}
		 AtomContainerManipulator.suppressHydrogens(molecule);
		return false;
	}

	/**
	 * This function will return the HashMap<aminoAcid_name, SMIRKS> for the 20 amino acids used in human body
	 * Reference: https://en.wikipedia.org/wiki/Amino_acid
	 * @return
	 * @throws Exception
	 */
	public HashMap<String, String> getAminoAcidMap(){
		HashMap<String, String> resultMap = new HashMap<>();
		//resultMap.put("Alanine", "[#6]-[#6](-[#7])-[#6](-[#8])=O");
		resultMap.put("Alanine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[H]");
		//resultMap.put("Arginine", "[#7]-[#6](-[#6]-[#6]-[#6]-[#7]-[#6](-[#7])=[#7])-[#6](-[#8])=O");
		resultMap.put("Arginine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])C([H])([H])C([H])([H])[#7]([H])-[#6](=[#7]\\[H])\\[#7]([H])[H]");
		//resultMap.put("Asparagine", "[#7]-[#6](-[#6]-[#6](-[#7])=O)-[#6](-[#8])=O");
		resultMap.put("Asparagine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#6](=O)-[#7]([H])[H]");
		//resultMap.put("Aspartate", "[#7]-[#6](-[#6]-[#6](-[#8])=O)-[#6](-[#8])=O");
		resultMap.put("Aspartate", "[H][#8]-[#6](=O)C([H])([H])C([H])([#7]([H])[H])[#6](=O)-[#8][H]");
		//resultMap.put("Cysteine", "[#7]-[#6](-[#6]-[#16])-[#6](-[#8])=O");
		resultMap.put("Cysteine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#16][H]");
		//resultMap.put("Glutamine", "[#7]-[#6](-[#6]-[#6]-[#6](-[#7])=O)-[#6](-[#8])=O");
		resultMap.put("Glutamine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])C([H])([H])[#6](=O)-[#7]([H])[H]");
		//resultMap.put("Glutamate", "[#7]-[#6](-[#6]-[#6]-[#6](-[#8])=O)-[#6](-[#8])=O");
		resultMap.put("Glutamate", "[H][#8]-[#6](=O)C([H])([H])C([H])([H])C([H])([#7]([H])[H])[#6](=O)-[#8][H]");
		//resultMap.put("Glycine", "[#7]-[#6]-[#6](-[#8])=O");
		resultMap.put("Glycine", "[H][#8]-[#6](=O)C([H])([H])[#7]([H])[H]");
		//resultMap.put("Histidine", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#7]-[#6]=[#7]-1)-[#6](-[#8])=O");
		resultMap.put("Histidine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#6]-1=[#6]([H])-[#7]([H])-[#6]([H])=[#7]-1");
		//resultMap.put("Isoleucine", "[#6]-[#6]-[#6](-[#6])-[#6](-[#7])-[#6](-[#8])=O");
		resultMap.put("Isoleucine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])(C([H])([H])[H])C([H])([H])C([H])([H])[H]");
		//resultMap.put("Leucine", "[#6]-[#6](-[#6])-[#6]-[#6](-[#7])-[#6](-[#8])=O");
		resultMap.put("Leucine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]");
		//resultMap.put("Lysine", "[#7]-[#6]-[#6]-[#6]-[#6]-[#6](-[#7])-[#6](-[#8])=O");
		resultMap.put("Lysine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[#7]([H])[H]");
		//resultMap.put("Methionine", "[#6]-[#16]-[#6]-[#6]-[#6](-[#7])-[#6](-[#8])=O");
		resultMap.put("Methionine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])C([H])([H])[#16]C([H])([H])[H]");
		//resultMap.put("Phenylalanine", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1)-[#6](-[#8])=O");
		resultMap.put("Phenylalanine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#6]-1=[#6]([H])-[#6]([H])=[#6]([H])-[#6]([H])=[#6]-1[H]");
		//resultMap.put("Proline", "[#8]-[#6](=O)-[#6]-1-[#6]-[#6]-[#6]-[#7]-1");
		resultMap.put("Proline", "[H][#8]-[#6](=O)C1([H])[#7]([H])C([H])([H])C([H])([H])C1([H])[H]");
		//resultMap.put("Serine", "[#7]-[#6](-[#6]-[#8])-[#6](-[#8])=O");
		resultMap.put("Serine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#8][H]");
		//resultMap.put("Theonine", "[#6]-[#6](-[#8])-[#6](-[#7])-[#6](-[#8])=O");		
		resultMap.put("Theonine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([#8][H])C([H])([H])[H]");	
		//resultMap.put("Tryptophan", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#7]-[#6]-2=[#6]-1-[#6]=[#6]-[#6]=[#6]-2)-[#6](-[#8])=O");
		resultMap.put("Tryptophan", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#6]-1=[#6]([H])-[#7]([H])-[#6]-2=[#6]-1-[#6]([H])=[#6]([H])-[#6]([H])=[#6]-2[H]");
		//resultMap.put("Tysosine", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#6]=[#6](-[#8])-[#6]=[#6]-1)-[#6](-[#8])=O");
		resultMap.put("Tysosine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])([H])[#6]-1=[#6]([H])-[#6]([H])=[#6](-[#8][H])-[#6]([H])=[#6]-1[H]");
		//resultMap.put("Valine", "[#6]-[#6](-[#6])-[#6](-[#7])-[#6](-[#8])=O");
		resultMap.put("Valine", "[H][#8]-[#6](=O)C([H])([#7]([H])[H])C([H])(C([H])([H])[H])C([H])([H])[H]");
		return resultMap;
	}
	

	
	public HashMap<String, String> getAminoAcidsForPeptide_old(){
		HashMap<String, String> resultMap = new HashMap<>();
		resultMap.put("Alanine", "[#6]-[#6](-[#7])-[#6]=O");
		resultMap.put("Arginine", "[#7]-[#6](-[#6]-[#6]-[#6]-[#7]-[#6](-[#7])=[#7])-[#6]=O");
		resultMap.put("Asparagine", "[#7]-[#6](-[#6]-[#6](-[#7])=O)-[#6]=O");		
		resultMap.put("Aspartate", "[#7]-[#6](-[#6]-[#6](-[#8])=O)-[#6]=O");	
		resultMap.put("Cysteine", "[#7]-[#6](-[#6]-[#16])-[#6]=O");		
		resultMap.put("Glutamine", "[#7]-[#6](-[#6]-[#6]-[#6](-[#7])=O)-[#6]=O");		
		resultMap.put("Glutamate", "[#7]-[#6](-[#6]-[#6]-[#6](-[#8])=O)-[#6]=O");
		resultMap.put("Glycine", "[#7]-[#6]-[#6]=O");		
		resultMap.put("Histidine", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#7]-[#6]=[#7]-1)-[#6]=O");		
		resultMap.put("Isoleucine", "[#6]-[#6]-[#6](-[#6])-[#6](-[#7])-[#6]=O");		
		resultMap.put("Leucine", "[#6]-[#6](-[#6])-[#6]-[#6](-[#7])-[#6]=O");		
		resultMap.put("Lysine", "[#7]-[#6]-[#6]-[#6]-[#6]-[#6](-[#7])-[#6]=O");		
		resultMap.put("Methionine", "[#6]-[#16]-[#6]-[#6]-[#6](-[#7])-[#6]=O");		
		resultMap.put("Phenylalanine", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1)-[#6]=O");		
		resultMap.put("Proline", "O=[#6]-[#6]-1-[#6]-[#6]-[#6]-[#7]-1");		
		resultMap.put("Serine", "[#7]-[#6](-[#6]-[#8])-[#6]=O");		
		resultMap.put("Theonine", "[#6]-[#6](-[#8])-[#6](-[#7])-[#6]=O");				
		resultMap.put("Tryptophan", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#7]-[#6]-2=[#6]-1-[#6]=[#6]-[#6]=[#6]-2)-[#6]=O");		
		resultMap.put("Tysosine", "[#7]-[#6](-[#6]-[#6]-1=[#6]-[#6]=[#6](-[#8])-[#6]=[#6]-1)-[#6]=O");		
		resultMap.put("Valine", "[#6]-[#6](-[#6])-[#6](-[#7])-[#6]=O");		
		return resultMap;
	}
	
	
	
}
