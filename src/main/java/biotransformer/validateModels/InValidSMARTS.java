package biotransformer.validateModels;

import java.util.ArrayList;
import java.util.HashMap;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.utils.HandlePolymers;

//public static void main(String[] args) throws Exception{
//	
//}

public class InValidSMARTS {	
	public HashMap<String, String> invalidSubstructuresMap = new HashMap<>();
	HandlePolymers hp = new HandlePolymers();
	public InValidSMARTS(){
		this.invalidSubstructuresMap = getInvalidStructureMap();
	}
	/**
	 * This function checks if the generated metabolite contains some obviously unreasonable substructure
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public boolean containInvalidSubstructure(IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		if(this.hp.isPolymer_moreThanDimer(oneMole)) return true;
		boolean contains = false;
		for(String key : this.invalidSubstructuresMap.keySet()){
			String smarts = this.invalidSubstructuresMap.get(key);
			Pattern pattern = SmartsPattern.create(smarts);
			if(pattern.matches(oneMole)){
				contains = true;
				break;
			}
			
		}
		AtomContainerManipulator.suppressHydrogens(oneMole);
		return contains;
	}
	

	public HashMap<String, String> getInvalidStructureMap(){
		HashMap<String, String> resultMap = new HashMap<>();
		resultMap.put("COOO", "[#8]-[#6](-[#8])-[#8]-*");//One carbon connecting to 3 oxygen is very unstable
		resultMap.put("C=OOO", "[#8]-[#6](=O)-[#8]-*");//One carbon connecting to 3 oxygen is very unstable
		resultMap.put("COOOO", "[#8]C([#8])([#8])[#8]-*");//one carbon connecting to 4 oxygens seems impossible. 
		resultMap.put("HydroxylEpOxide", "[#8]-[#6]-1-[#6]-[#8]-1");
		resultMap.put("DihydroxylOnOnecarbon", "[H][#8]-[#6]-[#8][H]");
		resultMap.put("C=C=O","[#6]=C=O");
		resultMap.put("O=C1CO1", "O=[#6]-1-[#6]-[#8]-1");
		resultMap.put("O-O", "[#8]-[#8]");
		resultMap.put("O=O", "O=O");
		return resultMap;
	}

}
