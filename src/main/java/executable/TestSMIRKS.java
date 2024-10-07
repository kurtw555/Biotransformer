package executable;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;

public class TestSMIRKS {
	//String smirks = "[H:18][#6:16]~1~[#6,#8,#7,#16:15]~[#6,#8,#7,#16:14]~[#6,#8,#7,#16:13]~[#6,#8,#7,#16:12]~2~[#6,#8,#7,#16;A;R2:17]~1~[#6,#8,#7,#16:1]~[#6,#8,#7,#16;A;R1:2]~[#6,#8,#7,#16;A;R2:3]~1-,=[#6,#8,#7,#16;A;R2:4]-,=3-,=[#6,#8,#7,#16:5]-,=[#6,#8,#7,#16:6]-,=[#6,#8,#7,#16;A:7]-,=[#6,#8,#7,#16:8]-,=3-,=[#6,#8,#7,#16:9]~[#6,#8,#7,#16;A:10]~[#6,#8,#7,#16:11]~2~1>>[H:18][#8]-[#6:16]~1~[#6,#8,#7,#16:15]~[#6,#8,#7,#16:14]~[#6,#8,#7,#16:13]~[#6,#8,#7,#16:12]~2~[#6,#8,#7,#16:11]~3~[#6,#8,#7,#16;A:10]~[#6,#8,#7,#16:9]-,=[#6,#8,#7,#16:8]-,=4-,=[#6,#8,#7,#16;A:7]-,=[#6,#8,#7,#16:6]-,=[#6,#8,#7,#16:5]-,=[#6,#8,#7,#16;A;R2:4]-,=4-,=[#6,#8,#7,#16;A;R2:3]~3~[#6,#8,#7,#16;A;R1:2]~[#6,#8,#7,#16:1]~[#6,#8,#7,#16;A;R2:17]~1~2";
	//String smirks = "[H:19][#8:5]-[#6:4]-1=[#6:15]-[#6:16]=[#6:17]-[#6:18]=[#6:3]-1-[#6:2](=[O:1])-[#6:7](=[O:8])-[#6:6]-[#6:9]-1=[#6:14]-[#6:13]=[#6:12]-[#6:11]=[#6:10]-1>>[H:19][#8:8][C:7]1([#6:6]-[#6:9]-2=[#6:14]-[#6:13]=[#6:12]-[#6:11]=[#6:10]-2)[#8:5]-[#6:4]-2=[#6:15]-[#6:16]=[#6:17]-[#6:18]=[#6:3]-2-[#6:2]1=[O:1]";
	String smirks = "[H:19][#8:5]-[c:4]1[c:15][c:16][c:17][c:18][c:3]1-[#6:2](=[O:1])-[#6:7](=[O:8])-[#6:6]-[c:9]1[c:10][c:11][c:12][c:13][c:14]1>>[H:19][#8:8][C:7]1([#6:6]-[c:9]2[c:10][c:11][c:12][c:13][c:14]2)[#8:5]-[c:4]2[c:15][c:16][c:17][c:18][c:3]2-[#6:2]1=[O:1]";
	//String smirks = "[H:19][#8:5]-[#6:4]=,:1=,:[#6:15]=,:[#6:16]=,:[#6:17]=,:[#6:18]=,:[#6:3]=,:1-[#6:2](=[O:1])-[#6:7](=[O:8])-[#6:6]-[#6:9]=,:1=,:[#6:10]=,:[#6:11]=,:[#6:12]=,:[#6:13]=,:[#6:14]=,:1>>[H:19][#8:8][C:7]1([#6:6]-[#6:9]-2=[#6:14]-[#6:13]=[#6:12]-[#6:11]=[#6:10]-2)[#8:5]-[#6:4]-2=[#6:15]-[#6:16]=[#6:17]-[#6:18]=[#6:3]-2-[#6:2]1=[O:1]";
	
	public static void main(String[] args) throws Exception{
//		String smiles = "CC(CCC(O)=O)C1CCC2C3CC(O)C4CC(O)CCC4(C)C3CC(O)C12C";
//		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//		IAtomContainer oneMole = sp.parseSmiles(smiles);
		TestSMIRKS ts = new TestSMIRKS();
//		ts.generateMetabolite(oneMole);
		//String smiles = "CCC(=O)N(C1CCN(CCC2=CC=CC=C2)CC1)C1=C(O)C=CC=C1"; //1.695893706516446
		//String smiles = "CCC(=O)N(C1CCN(CCC2=CC=CC=C2)CC1)C1=CC=CC=C1"; //1.841075533285448
		//String smiles = "Oc1cc(O)c(C(=O)C(=O)Cc2ccc(O)c(O)c2)c(O)c1"; //1.6348055121389788
		//String smiles = "OC1=CC(O)=C(C(=O)C(=O)CC2=CC=C(O)C(O)=C2)C(O)=C1";
		String smiles = "O[C@@H]1[C@@H](O)[C@@H](OC(O)[C@H]1O)CO";
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer oneMole = sp.parseSmiles(smiles);
		//PredictLogD logD_predictor = new PredictLogD();
		//System.out.println(logD_predictor.predictLogD(oneMole));
		//ts.generateMetabolite(oneMole);
		ts.checkMatch(oneMole);
	}
	
	public void generateMetabolite(IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(oneMole);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(oneMole);
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());		
		smrkMan.setFlagCheckAromaticityOnResultProcess(true);
		smrkMan.setFlagClearAromaticityBeforeResultProcess(true);
		SMIRKSReaction oneReaction = smrkMan.parse(this.smirks);
		
		IAtomContainerSet result = smrkMan.applyTransformationWithSingleCopyForEachPos(oneMole, null, oneReaction);
		for(int i = 0; i < result.getAtomContainerCount(); i++){
			System.out.println(sg.create(result.getAtomContainer(i)));
		}
	}
	public void checkMatch(IAtomContainer molecule) throws Exception{
		//String smarts = "[#8]-[#6]-[#6@@H]-1-[#8]-[#6](-[#8])-[#6@@H](-[#8])-[#6@H](-[#8])-[#6@H]-1-[#8]";
		//String smarts = "[#6]-1-[#6]-[#6]-[#8]-[#6]-[#6]-1";
		String smarts = "[#6]-1-[#6]-[#6]-[#6]-[#6]-[#6]-1";
		Pattern smp2 = SmartsPattern.create(smarts);
		boolean status = smp2.matches(molecule);
		System.out.println(status);
	}
}
