package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmilesGenerator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Biotransformer.bType;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;

public class BiotransformerSequence {
	//	public enum sequenceMode {
	//	
	//}
	
	protected ArrayList<BiotransformerSequenceStep> sequence;
	protected double scoreThreshold;
	protected Biotransformer utilityBiotransformer = null;
	boolean useDB;
	boolean useSubstitution;
	public BiotransformerSequence(ArrayList<BiotransformerSequenceStep> mySequence, boolean useDB, double scoreThreshold, boolean useSubstituition) {
		this.sequence 		= mySequence;
		this.scoreThreshold = scoreThreshold;
		this.useDB = useDB;
		this.useSubstitution = useSubstituition;
	}	
	public BiotransformerSequence(ArrayList<BiotransformerSequenceStep> mySequence, boolean useDB, boolean useSubstituition) {
		this.sequence = mySequence;
		this.scoreThreshold = 0.0;
		this.useDB = useDB;
		this.useSubstitution = useSubstituition;
	}

	
	public BiotransformerSequence(String mySequence, boolean useDB, double scoreThreshold, boolean useSubstituition) throws Exception {
		this.sequence 		= createSequenceFromString(mySequence, scoreThreshold);
		this.scoreThreshold = scoreThreshold;
		this.useDB = useDB;
		this.useSubstitution = useSubstituition;
	}
	public BiotransformerSequence(String mySequence, boolean useDB, boolean useSubstituition) throws Exception {
		this.sequence 		= createSequenceFromString(mySequence, 0.0);
		this.scoreThreshold = 0.0;
		this.useDB = useDB;
		this.useSubstitution = useSubstituition;
	}	
	
	
	
	private ArrayList<BiotransformerSequenceStep> createSequenceFromString(String mySequence, double scoreThreshold) throws Exception{
		ArrayList<BiotransformerSequenceStep> final_seq = new ArrayList<BiotransformerSequenceStep>();
		try {
			ArrayList<BiotransformerSequenceStep> seq = new ArrayList<BiotransformerSequenceStep>();
			String[] steps = mySequence.split(";");
			for(String step : steps) {
				String[] attr = step.split(":");
				if(attr.length==2) {
					String btype_str = attr[0].trim().toUpperCase();
					int nsteps = Integer.valueOf(attr[1].trim());
					
					try {
						bType btype_ = Biotransformer.bType.valueOf(btype_str);
					}
					catch (IllegalArgumentException ie) {
						throw new IllegalArgumentException("Illegal biotransformer sequence type '" + btype_str + "'.");
					}
					catch (Exception e) {
						System.err.println(e.getMessage());
					}							
					seq.add(new BiotransformerSequenceStep(Biotransformer.bType.valueOf(btype_str), nsteps, scoreThreshold));
				}
				else {
					throw new IllegalArgumentException("Illegal biotransformer sequence '" + mySequence + "\n\tMake sure your sequence is semi-colon separated, and that each step is colon-separated pair of one valid biotransformer type, and one integer.");					
					
				}
			}	
			final_seq = seq;
		}
		catch (Exception e){
			System.out.println(e);
		}

		return final_seq;
	}
	
	public String toString() {
		String str_representation = "";
		ArrayList<String> steps = new ArrayList<String>();
		int counter = 0;
		for(BiotransformerSequenceStep s : this.sequence) {
			counter++;
			steps.add("Step " + counter + ": " + s.toString());
		}
		str_representation = String.join("\n", steps);
		return str_representation;
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainer startingCompound, double scoreThreshold) throws Exception{
		return runSequence(startingCompound, scoreThreshold, 1);
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainer startingCompound, double scoreThreshold, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		LinkedHashMap<Biotransformer.bType, Object> btransformers = new LinkedHashMap<Biotransformer.bType, Object>();
		IAtomContainerSet currentMetabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		currentMetabolites.addAtomContainer(startingCompound);
		
		for(BiotransformerSequenceStep step : this.sequence) {
			ArrayList<Biotransformation> currentBiots =  new ArrayList<Biotransformation>();
				
			if(step.btype == bType.CYP450) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new Cyp450BTransformer(BioSystemName.HUMAN, this.useDB, this.useSubstitution));
				}
				currentBiots = ((Cyp450BTransformer) btransformers.get(step.btype)).predictCyp450BiotransformationChainByMode(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold, cyp450Mode);			
			}
			else if(step.btype == bType.ECBASED) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new ECBasedBTransformer(BioSystemName.HUMAN, this.useDB, this.useSubstitution));
				}
				//currentBiots = ((ECBasedBTransformer) btransformers.get(step.btype)).simulateECBasedMetabolismChain(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);		
				currentBiots = ((ECBasedBTransformer) btransformers.get(step.btype)).simulateECBasedMetabolismChain(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);	
			}
			else if(step.btype == bType.ENV) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new EnvMicroBTransformer());
				}
				currentBiots = ((EnvMicroBTransformer) btransformers.get(step.btype)).applyEnvMicrobialTransformationsChain(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);
			}			
			else if(step.btype == bType.HGUT) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new HGutBTransformer(this.useDB, this.useSubstitution));
				}
				currentBiots = ((HGutBTransformer) btransformers.get(step.btype)).simulateGutMicrobialMetabolism_withDepolymerization(currentMetabolites, true, true, step.nOfIterations, this.scoreThreshold);						
			}
			else if(step.btype == bType.PHASEII) {
				if(btransformers.get(step.btype) == null) {
					btransformers.put(step.btype, new Phase2BTransformer(BioSystemName.HUMAN, this.useDB, this.useSubstitution));
				}
				currentBiots = ((Phase2BTransformer) btransformers.get(step.btype)).applyPhase2TransformationsChainAndReturnBiotransformations(currentMetabolites, true, true, true, step.nOfIterations, this.scoreThreshold);										
			}
			
			currentMetabolites.add(((Biotransformer) btransformers.get(step.btype)).extractProductsFromBiotransformations(currentBiots));
			biotransformations.addAll(currentBiots);
		}
		
//		this.utilityBiotransformer = ((Collection<Entry<bType, Object>>) btransformers.entrySet()).stream().reduce((first, second) -> second).orElse(null)biotransformations.gtValue();
//		for(Biotransformation b: biotransformations) {
//			b.display();		
//		}
		

		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	

	public ArrayList<Biotransformation> runSequence(IAtomContainerSet startingCompounds, double scoreThreshold) throws Exception{
		return runSequence(startingCompounds, scoreThreshold, 1);
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainerSet startingCompounds, double scoreThreshold, int cyp450Mode) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer starting_ac : startingCompounds.atomContainers()) {
			try {
				biotransformations.addAll(runSequence(starting_ac, scoreThreshold, cyp450Mode));
			}
			catch (Exception e) {
				System.out.println(e);
			}
			
		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
}
