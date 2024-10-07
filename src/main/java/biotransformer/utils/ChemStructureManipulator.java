/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.util.ArrayList;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
//import ambit2.tautomers.*;
//import ambit2.tautomers.processor.StructureStandardizer;
import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MetabolicReaction;

public class ChemStructureManipulator {
	protected static SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//	protected static StructureStandardizer sstandardizer = new ambit2.tautomers.processor.StructureStandardizer();
	
	public ChemStructureManipulator() {
		// TODO Auto-generated constructor stub
		smrkMan.setFlagApplyStereoTransformation(false);
		smrkMan.setFlagCheckResultStereo(true);
		smrkMan.setFlagFilterEquivalentMappings(true);
		smrkMan.setFlagProcessResultStructures(true);
		smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);
		
//		sstandardizer.setEnabled(true);
//		sstandardizer.setImplicitHydrogens(true);
//		sstandardizer.setGenerate2D(true);
//		sstandardizer.setNeutralise(false);
////		sstandardizer.setGenerateStereofrom2D(true);
//		sstandardizer.setGenerateSMILES(true);	
//		sstandardizer.setGenerateTautomers(true);
	}
	
	/**
	 * This function applies some preprocessing operations, such as setting the
	 * flag of atoms from aromatic rings to "ISAROMATIC", and kelulizing
	 * molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return : A processed molecule (AtomContainer)
	 * @throws Exception 
	 */
//	public static IAtomContainer preprocessContainer(IAtomContainer molecule)
//			throws Exception {
//
//		IAtomContainer molClone =  molecule.clone();
//
//		return sstandardizer.process(molClone);	
//	}
	
	public static IAtomContainer preprocessContainer(IAtomContainer molecule)
			throws CDKException {
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
//	    Aromaticity aromaticity = new Aromaticity( ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));
		
		for (IBond bond : molecule.bonds()) {
			if (bond.isAromatic() && bond.getOrder() == IBond.Order.UNSET) {
				bond.setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);

			} 
		}
		aromaticity.apply(molecule);
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(molecule);
		sdg.generateCoordinates();		
		IAtomContainer layedOutMol = sdg.getMolecule();

		return layedOutMol;
	}
	public static IAtomContainer standardizeMoleculeWithCopy(IAtomContainer molecule) throws Exception{
		return  standardizeMoleculeWithCopy(molecule, true);
	}

	public static IAtomContainer standardizeMoleculeWithCopy(IAtomContainer molecule, boolean preprocess) throws Exception{
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		for(int k = 0; k < molecule.getAtomCount(); k++) {
			IAtom oneAtom = molecule.getAtom(k);
			oneAtom.setProperty("AtomIdx", k);
		}
		IAtomContainer molClone =  molecule.clone();
		if(preprocess){
			molClone = ChemStructureManipulator.preprocessContainer(molClone);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(molClone);
		} else {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molClone);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(molClone);
		}
		
		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
		MReactionSets mrs = new MReactionSets();
		//for(MetabolicReaction mr : MReactionSets.standardizationReactions){
		for(MetabolicReaction mr : mrs.standardizationReactions){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(mr, molClone)){
				matchedReactions.add(mr);
			}
		}				
		for(MetabolicReaction m : matchedReactions){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, molClone)){
				smrkMan.applyTransformation(molClone, m.getSmirksReaction());
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molClone);
				SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
				String smiles = sg.create(molClone);
				if(smiles.contains("*")) {
					molClone = getCorrectMetabolite(molecule,molClone);
				}		
				adder.addImplicitHydrogens(molClone);			
			}

		}

		AtomContainerManipulator.suppressHydrogens(molClone);
		return molClone;
	}


	public static IAtomContainer getCorrectMetabolite(IAtomContainer reactant, IAtomContainer metabolite) throws Exception{
		for(int k = 0; k < metabolite.getAtomCount(); k++) {
			if(metabolite.getAtom(k).getSymbol().equalsIgnoreCase("R") || metabolite.getAtom(k).getSymbol().equalsIgnoreCase("A")) {
				Integer checkIdx = metabolite.getAtom(k).getProperty("AtomIdx");
				for(int j = 0; j < reactant.getAtomCount(); j++) {
					if(reactant.getAtom(j).getProperty("AtomIdx") != null && reactant.getAtom(j).getProperty("AtomIdx") == checkIdx) {
						metabolite.getAtom(k).setAtomicNumber(reactant.getAtom(j).getAtomicNumber());
						break;
					}
				}
			}
		}
		return metabolite;
	}
	
	

}
