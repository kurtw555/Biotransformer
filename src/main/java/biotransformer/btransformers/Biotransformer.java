
/**
 * This class implements several functions that apply transformations to
 * molecules, collect resulting products, and store them in various data
 * structures or formats. I uses several other classes such as Enzymes, and
 * MetabolicReactions.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */


package biotransformer.btransformers;



import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import org.apache.commons.lang3.StringUtils;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.JsonParser.Feature;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.dbrelevant.RetriveFromDB;
import biotransformer.esaprediction.ESSpecificityPredictor;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MReactionsFilter;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.ChemdbRest;
import biotransformer.utils.ChemicalClassFinder;
import biotransformer.utils.ChemicalClassFinder.ChemicalClassName;
import exception.BioTransformerException;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.HandlePolymers;
import biotransformer.utils.Utilities;
import biotransformer.validateModels.InValidSMARTS;




/**
 * @author Djoumbou Feunang, Yannick
 *
 */

public class Biotransformer {

	public enum bType {
		ALLHUMAN, CYP450, ECBASED, ENV, HGUT, PHASEII, SUPERBIO, ABIOTIC
	}
	public boolean useDB;
	public boolean useSubstitution;
	public RetriveFromDB rfdb = null;
	InValidSMARTS invalidSMARTS = new InValidSMARTS();
	protected SMIRKSManager smrkMan;
	public BioSystem bSystem;
	protected MReactionsFilter mRFilter;
	protected LinkedHashMap<String, Double> reactionORatios;
	protected IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
	public SmilesGenerator smiGen 		= new SmilesGenerator(SmiFlavor.Isomeric);
	public boolean inAllHuman = false;
	
	protected SmilesParser	smiParser		= new SmilesParser(builder);
	public InChIGeneratorFactory inchiGenFactory;
	public LinkedHashMap<String, ArrayList<MetabolicReaction>>	reactionsByGroups 
											= new LinkedHashMap<String, ArrayList<MetabolicReaction>>();
	public LinkedHashMap<String, ArrayList<Enzyme>>	enzymesByreactionGroups 
											= new LinkedHashMap<String, ArrayList<Enzyme>>();
		
	public ArrayList<Enzyme> enzymesList = new ArrayList<Enzyme>();
	public LinkedHashMap<String, MetabolicReaction> reactionsHash = new LinkedHashMap<String, MetabolicReaction>();
	
	public ObjectMapper mapper = new ObjectMapper();
	protected ESSpecificityPredictor esspredictor;
	public HandlePolymers hp = new HandlePolymers();
		
	public Biotransformer(BioSystemName bioSName, boolean useDB, boolean useSubstitution) throws JsonParseException, JsonMappingException, FileNotFoundException, IOException, BioTransformerException, CDKException{
		this.useDB = useDB;
		this.useSubstitution = useSubstitution;
		this.mapper.configure(Feature.ALLOW_COMMENTS, true);
		this.mapper.configure(Feature.ALLOW_BACKSLASH_ESCAPING_ANY_CHARACTER, true);
		this.bSystem = new BioSystem(bioSName, mapper);		
		this.smrkMan = this.bSystem.getSmirksManager();
		this.esspredictor = new ESSpecificityPredictor(this.bSystem);
		this.inchiGenFactory = InChIGeneratorFactory.getInstance();
		
		this.smrkMan.setFlagApplyStereoTransformation(false);
		this.smrkMan.setFlagCheckResultStereo(true);
		this.smrkMan.setFlagFilterEquivalentMappings(true);
		this.smrkMan.setFlagProcessResultStructures(true);
		this.smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);

				
		this.reactionORatios = this.bSystem.getReactionsORatios();
		setReactionsGroups();
		setReactionsFilter(); 
		
		
	}
	/**
	 * This function will explore all metabolic pathways for the query compound stored in HMDB and return the results
	 * @param substrates
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> getBiotransoformationsByExploringHDMB(IAtomContainer substrate) throws Exception{
		ArrayList<String> substrateProcessed = new ArrayList<>();	
		ArrayList<Biotransformation> results = new ArrayList<>();
		ArrayList<Biotransformation> currentBiotransformations = this.rfdb.getBiotransformationRetrievedFromDB(substrate, true);
		results.addAll(currentBiotransformations);
		IAtomContainerSet substratePool = extractProductsFromBiotransformations(currentBiotransformations);	
		//substratePool = Utilities.getNovelSubstrates(substratePool, substrateProcessed);
		while(substratePool != null && !substratePool.isEmpty()) {
			currentBiotransformations = new ArrayList<>();
			for(int i = 0; i < substratePool.getAtomContainerCount(); i++) {
				IAtomContainer oneMole = substratePool.getAtomContainer(i);
				if(countNonHydrogenAtoms(oneMole) <= 4) continue;
				ArrayList<Biotransformation> db_results = this.rfdb.getBiotransformationRetrievedFromDB(oneMole, true);
				for(int t = 0; t < db_results.size(); t++) {
					if(isConjugate(db_results.get(t))) {
						db_results.remove(t);
						t = t-1;
					}
				}
				currentBiotransformations.addAll(db_results);				
			}
			results.addAll(currentBiotransformations);
			substratePool = extractProductsFromBiotransformations(currentBiotransformations);	
			substratePool = Utilities.getNovelSubstrates(substratePool, substrateProcessed);		
			substrateProcessed = Utilities.updateProcessedSubstratePool(substrateProcessed, substratePool);
		}
		return Utilities.selectUniqueBiotransformations(results);
		//return Utilities.selectByUniqueMetabolites(results);
	}
	public boolean isConjugate(Biotransformation bt) throws Exception{		
		IAtomContainer substrate = bt.getSubstrates().getAtomContainer(0);
		IAtomContainerSet metabolites = bt.getProducts();
		IAtomContainer metabolite = metabolites.getAtomContainer(0);
		if(metabolite.getProperty("Reaction ID").equals(2365)) {

		}
		int currentCounter = getNonHydrogenAtoms(metabolite);
		for(int i = 0; i < metabolites.getAtomContainerCount(); i++) {
			if(getNonHydrogenAtoms(metabolites.getAtomContainer(i)) > currentCounter) {
				currentCounter = getNonHydrogenAtoms(metabolites.getAtomContainer(i));
				metabolite = metabolites.getAtomContainer(i);
			}
		}
		if((currentCounter - getNonHydrogenAtoms(substrate)) > 3) return true;
		else return false;
		
	}
	public int getNonHydrogenAtoms(IAtomContainer metabolite) {
		int counter = 0;
		for(int t = 0; t < metabolite.getAtomCount(); t++) {
			if(!metabolite.getAtom(t).getSymbol().equalsIgnoreCase("H")) counter++;
		}
		return counter;
	}
	private void setReactionsGroups() throws JsonParseException, JsonMappingException, FileNotFoundException, IOException{
		MReactionSets mrs = new MReactionSets();
		//reactionsByGroups.put("standardizationReactions", MReactionSets.standardizationReactions);
		reactionsByGroups.put("standardizationReactions", mrs.standardizationReactions);
		MReactionSets msets = new MReactionSets();
		for(MetabolicReaction mr : msets.standardizationReactions) {
			this.reactionsHash.put(mr.name, mr);
		}
	}
	
	private void setReactionsFilter(){
		this.mRFilter = new MReactionsFilter(this.bSystem);
	}
	
	public BioSystemName getBioSystemName(){
		return this.bSystem.name;
	}
	
	
	public LinkedHashMap<String, ArrayList<MetabolicReaction>> getReactionsList(){
		return this.reactionsByGroups;
	}
	


	/**
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @param mReaction
	 *            : A MetabolicReaction
	 * @param preprocess
	 *            : A boolean that specifies whether to apply preprocessing or
	 *            not
	 * @return : A set of metabolites or products generated by applying the
	 *         reaction to the molecule of interest
	 * @throws Exception
	 *  		  : Throws an Exception
	 */
	public IAtomContainerSet generateAllMetabolitesFromAtomContainer(IAtomContainer molecule, MetabolicReaction mReaction,boolean preprocess) throws Exception {
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		metabolites = generateAllMetabolitesFromAtomContainer(molecule, mReaction.getSmirksReaction(), preprocess);	
		return metabolites;
	}

	/**
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @param reaction
	 *            : A SMIRKSReaction
	 * @param preprocess
	 *            : A boolean that specifies whether to apply preprocessing or
	 *            not
	 * @return : A set of metabolites or products generated by applying the
	 *         reaction to the molecule of interest
	 * @throws Exception
	 *  			  : Throws an Exception
	 */
	public IAtomContainerSet generateAllMetabolitesFromAtomContainer(IAtomContainer molecule, SMIRKSReaction reaction,boolean preprocess) throws Exception {
		// https://github.com/ideaconsult/examples-ambit/tree/master/smirks-example

		IAtomContainer reactant = (IAtomContainer) molecule.clone();
		if (preprocess) {
			try {
				reactant = ChemStructureManipulator.preprocessContainer(reactant);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(reactant);
			}
			catch (Exception e) {
				System.out.println(e);
			}
		} 
		else {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(reactant);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(reactant);
		}
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(reactant);
		for(int k = 0; k < reactant.getAtomCount(); k++) {
			IAtom oneAtom = reactant.getAtom(k);
			oneAtom.setProperty("AtomIdx", k);
		}
		IAtomContainerSet metabolites = this.smrkMan.applyTransformationWithSingleCopyForEachPos(reactant, null, reaction);
		IAtomContainerSet postprocessed_metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		if (metabolites != null) {
			int nr_of_metabolites = metabolites.getAtomContainerCount();
			if (nr_of_metabolites > 0) {
				for (int i = 0; i < nr_of_metabolites; i++) {
					IAtomContainerSet partitions = ChemStructureExplorer.checkConnectivity(metabolites.getAtomContainer(i));
					for (IAtomContainer c :  partitions.atomContainers()) {
						if(!isValidMetabolte(c)) continue;
						if(!ChemStructureExplorer.isUnneccessaryMetabolite(c)){
//							aromaticity.apply(c);
							try{
								postprocessed_metabolites.addAtomContainer(ChemStructureManipulator.preprocessContainer(c));
							}catch (NullPointerException n) {
								System.err.println(n.getMessage());
								postprocessed_metabolites.addAtomContainer(c);
							}
							
						}
					}
				}
			}
		}
		IAtomContainerSet final_metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < postprocessed_metabolites.getAtomContainerCount(); i++) {
			IAtomContainer oneMetabolite = postprocessed_metabolites.getAtomContainer(i);
			String smiles = this.smiGen.create(oneMetabolite);
			if(smiles.contains("*")) {
				IAtomContainer corrected = ChemStructureManipulator.getCorrectMetabolite(reactant, oneMetabolite);
				corrected.addProperties(oneMetabolite.getProperties());
				final_metabolites.addAtomContainer(corrected);
			}
			else final_metabolites.addAtomContainer(oneMetabolite);
		}
		//return ChemStructureExplorer.uniquefy(postprocessed_metabolites);
		return ChemStructureExplorer.uniquefy(final_metabolites);
	}
	
	
	public IAtomContainerSet generateAllMetabolitesFromAtomContainerViaTransformationAtAllLocations(
			IAtomContainer molecule, SMIRKSReaction reaction,
			boolean preprocess) throws Exception {
		// https://github.com/ideaconsult/examples-ambit/tree/master/smirks-example

		IAtomContainer reactant = (IAtomContainer) molecule.clone();
		if (preprocess) {
			try {
				reactant = ChemStructureManipulator.preprocessContainer(reactant);
			}
			catch (Exception e) {
				System.out.println(e);
			}
		} 
		else {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(reactant);
//			AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		}
		IAtomContainerSet postprocessed_metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);

//		Aromaticity aromaticity = new Aromaticity( ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));
		if (reactant != null) {
			IAtomContainerSet partitions = ChemStructureExplorer.checkConnectivity(reactant);		
			
			for (int k = 0; k < partitions.getAtomContainerCount(); k++) {
				aromaticity.apply(partitions.getAtomContainer(k));
				postprocessed_metabolites.addAtomContainer(partitions.getAtomContainer(k));

			}
			
		}
		return ChemStructureExplorer.uniquefy(postprocessed_metabolites);
	}

	public ArrayList<Biotransformation> applyReactionsAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter) throws Exception{
		return applyReactionsAndReturnBiotransformations(target, reactions,preprocess, filter, 0.0);
	}
	
	public ArrayList<Biotransformation> applyReactionsAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();		
		IAtomContainer starget = ChemStructureManipulator.standardizeMoleculeWithCopy(target, preprocess);
	
		if(target.getProperty("InChI") == null){
			try {
			InChIGenerator gen0 = this.inchiGenFactory.getInChIGenerator(target);
			target.setProperty("InChI", gen0.getInchi());
			target.setProperty("InChIKey", gen0.getInchiKey());
			} catch (CDKException c){
				System.err.println(c.getLocalizedMessage());
			}	
		}

		Utilities.addPhysicoChemicalProperties(target);
		
		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
		ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();
		for (MetabolicReaction i : reactions) {		
//			Work on this. It caused errors for LacCer(d20-1_16-1(9Z)) and other glycosyceramides
			boolean match_constraints = ChemStructureExplorer.compoundMatchesReactionConstraints(i, starget);
			if (match_constraints) {				
				matchedReactions.add(i);
			}
		}		
		if(filter == false){
			filteredReactions = matchedReactions;		
		} else{
			filteredReactions = new ArrayList<MetabolicReaction>(this.mRFilter.filterReactions(matchedReactions).values());
		}

		for(MetabolicReaction j : filteredReactions){
			IAtomContainerSet partialSet = this.generateAllMetabolitesFromAtomContainer(starget, j, true);
			Double score=0.0;
			AtomContainerSet subs = new AtomContainerSet();
			AtomContainerSet prod = new AtomContainerSet();
					
			if(partialSet.getAtomContainerCount()>0){				
				if(target.getProperty("Score") != null){					
					score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(j.name)  );
				}else{
					score = this.bSystem.getReactionsORatios().get(j.name);
				}
				if(score>=scoreThreshold){
					subs.addAtomContainer(target);
					for(IAtomContainer pc : partialSet.atomContainers()){
						if(!ChemStructureExplorer.isUnneccessaryMetabolite(pc)){
//						AtomContainerManipulator.suppressHydrogens(pc);
							try{
							InChIGenerator gen = this.inchiGenFactory.getInChIGenerator(pc);
							pc.setProperty("InChI", gen.getInchi());
							pc.setProperty("InChIKey", gen.getInchiKey());
							}catch (CDKException c){
								System.err.println(c.getLocalizedMessage());
							}
							Utilities.addPhysicoChemicalProperties(pc);
							prod.addAtomContainer(AtomContainerManipulator.removeHydrogens(pc));
						}
					}
									
					ArrayList<Enzyme> enzList = new ArrayList<Enzyme>();
//					j.display();
					
					Biotransformation bioT = new Biotransformation(subs, j.name, null, prod, score, this.getBioSystemName());
					results.add(bioT);
				}
			}	

		}
		return results;
	}

	public ArrayList<Biotransformation> applyReactionAndReturnBiotransformations(IAtomContainer target,
			MetabolicReaction reaction, boolean preprocess) throws Exception{
		return applyReactionAndReturnBiotransformations(target, reaction, preprocess, 0.0);
	}	
	
	public ArrayList<Biotransformation> applyReactionAndReturnBiotransformations(IAtomContainer target,
			MetabolicReaction reaction, boolean preprocess, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();		
		IAtomContainer starget = ChemStructureManipulator.standardizeMoleculeWithCopy(target, preprocess);
	
		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
		
		InChIGenerator gen0 = this.inchiGenFactory.getInChIGenerator(target);
		target.setProperty("InChI", gen0.getInchi());
		target.setProperty("InChIKey", gen0.getInchiKey());
		target.setProperty("SMILES", this.smiGen.create(target));
		Utilities.addPhysicoChemicalProperties(target);
		target.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(target));
	
		boolean match_constraints = ChemStructureExplorer.compoundMatchesReactionConstraints(reaction, starget);
		if (match_constraints) {
			matchedReactions.add(reaction);
		}


		IAtomContainerSet partialSet = generateAllMetabolitesFromAtomContainer(starget, reaction, false);
		Double score=0.0;
		AtomContainerSet subs = new AtomContainerSet();
		AtomContainerSet prod = new AtomContainerSet();
		
		
		if(partialSet.getAtomContainerCount()>0){
			
			if(target.getProperty("Score") !=null){	
				
				score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(reaction.name)  );
			}else{
				score = this.bSystem.getReactionsORatios().get(reaction.name);
			}
	
			if(score>=scoreThreshold){
				subs.addAtomContainer(target);
				for(IAtomContainer pc : partialSet.atomContainers()){
//					AtomContainerManipulator.suppressHydrogens(pc);
					InChIGenerator gen = this.inchiGenFactory.getInChIGenerator(pc);
					pc.setProperty("InChI", gen.getInchi());
					pc.setProperty("InChIKey", gen.getInchiKey());
					pc.setProperty("SMILES", this.smiGen.create(pc));
					pc.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(pc));
					Utilities.addPhysicoChemicalProperties(pc);
					prod.addAtomContainer(AtomContainerManipulator.removeHydrogens(pc));
				}
				
				Biotransformation bioT = new Biotransformation(subs, reaction.name, null, prod, score, this.getBioSystemName() );
				results.add(bioT);
			}
		}	

		return results;
	}

	public ArrayList<Biotransformation> applyReactionsFromContainersAndReturnBiotransformations(IAtomContainerSet targets,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> products = new ArrayList<Biotransformation>();
		
		for(IAtomContainer ac : targets.atomContainers()){
			products.addAll(applyReactionsAndReturnBiotransformations(ac, reactions, preprocess, filter, scoreThreshold)); 
		}		
			
		return products;
	}
	
	public ArrayList<Biotransformation> applyReactionsChainAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		return applyReactionsChainAndReturnBiotransformations(target, reactions, preprocess, filter, 
				nr_of_steps, 0.0);
	}
	
	public ArrayList<Biotransformation> applyReactionsChainAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{

		AtomContainerSet startingSet = new AtomContainerSet();
		startingSet.addAtomContainer(target);
		
		return applyReactionsChainAndReturnBiotransformations(startingSet, reactions, 
				preprocess, filter, nr_of_steps, scoreThreshold);
	}
		
	public ArrayList<Biotransformation> applyReactionsChainAndReturnBiotransformations(IAtomContainerSet targets,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		return applyReactionsChainAndReturnBiotransformations(targets, reactions, preprocess, 
				filter, nr_of_steps, 0.0);
	}
	
	public ArrayList<Biotransformation> applyReactionsChainAndReturnBiotransformations(IAtomContainerSet targets,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> products = new ArrayList<Biotransformation>();		
//		IAtomContainerSet containers = (AtomContainerSet) targets;
		IAtomContainerSet containers = targets;
		
		int counter = 0;
		
		while(nr_of_steps>0){
			counter++;			
			ArrayList<Biotransformation> currentProducts = applyReactionsFromContainersAndReturnBiotransformations(containers, reactions, preprocess, filter, scoreThreshold);
			nr_of_steps--;
//			System.err.println(currentProducts.size() + " biotransformations at step " + counter);
			if(!currentProducts.isEmpty()){
				
				
				products.addAll(currentProducts);
				containers.removeAllAtomContainers();
				containers = extractProductsFromBiotransformations(currentProducts);
			}
			else{
				break;
			}
		}
		return products;
	}
	
	public IAtomContainerSet applyReactions(IAtomContainer target, SMIRKSManager smrkMan,
			LinkedHashMap<String, MetabolicReaction> reactions, boolean preprocess, boolean filter) throws Exception {
		return applyReactions(target,
				reactions, preprocess, filter, 0.0);
	}
	
	public IAtomContainerSet applyReactions(IAtomContainer target,
			LinkedHashMap<String, MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception {
		IAtomContainerSet products = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		int count = 0;
		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
		ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();
		IAtomContainer starget = target.clone();
		
		if (preprocess) {
			try {
//				starget = standardizeMoleculeWithCopy(target);
				starget = ChemStructureManipulator.preprocessContainer(target);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
			}
			catch (Exception e) {
				System.out.println(e);
			}
		}
		else{
//			starget = standardizeMoleculeWithCopy(target, false);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
			
		}
		
		InChIGenerator gen0 = this.inchiGenFactory.getInChIGenerator(target);
		

		for (MetabolicReaction i : reactions.values()) {
	
			boolean match_constraints = ChemStructureExplorer.compoundMatchesReactionConstraints(i, starget);
			
			if (match_constraints) {
				matchedReactions.add(i);
			}
		}
		
		
		if(filter == false){
			filteredReactions = matchedReactions;	
		} else{
			filteredReactions = new ArrayList<MetabolicReaction>(this.mRFilter.filterReactions(matchedReactions).values());
		}
		
		
		for(MetabolicReaction j : filteredReactions){
			IAtomContainerSet partialSet = generateAllMetabolitesFromAtomContainer(starget, j, false);
			if(partialSet.getAtomContainerCount()>0){
				for(IAtomContainer pc : partialSet.atomContainers()){
					
					Double score=0.0;
					
					if(target.getProperty("Score") !=null){		
						score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(j.name)  );
					}else{
						score = this.bSystem.getReactionsORatios().get(j.name);
					}
					
					if(score>=scoreThreshold){
						pc.setProperty("Reaction", j.name);
						InChIGenerator gen = this.inchiGenFactory.getInChIGenerator(pc);
						pc.setProperty("Precursor", gen0.getInchiKey());
						pc.setProperty("InChI", gen.getInchi());
						pc.setProperty("InChIKey", gen.getInchiKey());
						pc.setProperty("Score", score);
						Utilities.addPhysicoChemicalProperties(pc);
						products .addAtomContainer(AtomContainerManipulator.removeHydrogens(pc));
					}
				}
			}		
		}

		IAtomContainerSet unique = ChemStructureExplorer.uniquefy(products);
		
		return unique;
	}

	public IAtomContainerSet applyReactionChain(IAtomContainer target, SMIRKSManager smrkMan,
			LinkedHashMap<String, MetabolicReaction> reactions, boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		return applyReactionChain(target, smrkMan, reactions, preprocess, filter, nr_of_steps, 0.0);
	}
	
	public IAtomContainerSet applyReactionChain(IAtomContainer target, SMIRKSManager smrkMan,
			LinkedHashMap<String, MetabolicReaction> reactions, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		IAtomContainerSet products = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);


	    int i = nr_of_steps;	
		int step = 0;
		boolean stop = false;
		
		if (preprocess) {
			try {
				target = ChemStructureManipulator.preprocessContainer(target);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
			}
			catch (Exception e) {
				System.out.println(e);
			}
		}else{
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(target);
		}
		
		while (nr_of_steps > 0){
			step++;	
			IAtomContainerSet currentProducts = applyReactions(target, this.smrkMan, reactions, false, filter);
			nr_of_steps--;
			
			if(currentProducts.getAtomContainerCount() > 0){
				int j = 1;
				for( IAtomContainer ac : currentProducts.atomContainers() ){
					if(target.getProperty(CDKConstants.TITLE) == null){
						ac.setProperty(CDKConstants.TITLE,"metabolite_"+j);
					} else{
						ac.setProperty(CDKConstants.TITLE,target.getProperty(CDKConstants.TITLE)+"_"+j);
					}
					j++;

					products.addAtomContainer(ac);
					IAtomContainerSet t = applyReactionChain(ac, this.smrkMan, reactions, false, filter, nr_of_steps, scoreThreshold);
					
					if(t.getAtomContainerCount() > 0){
						products.add( t );
					}				
				}
			}
			else{
				stop = true;
				break;
			}
		}	
		System.err.println(products.getAtomContainerCount() + " products after " + step + " step(s).");			
		return ChemStructureExplorer.uniquefy(products);

	}	

	public ArrayList<Biotransformation> metabolizeWithEnzyme(IAtomContainer substrate, String enz, boolean preprocess, boolean filter, double threshold) throws Exception {
		return metabolizeWithEnzyme(substrate, enz, true, preprocess,filter, threshold);
	}
	
	public ArrayList<Biotransformation> metabolizeWithEnzyme(IAtomContainer substrate, String enz, boolean predictESSpecificity, boolean preprocess, boolean filter, double threshold) throws Exception {
		IAtomContainer clonedSubs = substrate.clone();
		
//		if(preprocess){
//			clonedSubs = ChemStructureManipulator.preprocessContainer(clonedSubs);
//			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clonedSubs);
//			AtomContainerManipulator.convertImplicitToExplicitHydrogens(clonedSubs);	
//		} else{
//			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clonedSubs);
//			AtomContainerManipulator.convertImplicitToExplicitHydrogens(clonedSubs);
//		}
		
		try{
			if(this.bSystem.getEnzymeHash().containsKey(enz)){
				Enzyme e = this.bSystem.getEnzymeHash().get(enz);		
				return metabolizeWithEnzyme(substrate, e, predictESSpecificity, preprocess, filter, threshold);
				
			} else {
				throw new IllegalArgumentException(enz.toString() + " is not associated with the biosystem " + this.getBioSystemName());
			}
		} catch (IllegalArgumentException iae) {
			System.err.println(iae.getLocalizedMessage());
			return null;
		}
	}

	public ArrayList<Biotransformation> metabolizeWithEnzyme(IAtomContainer substrate, Enzyme enzyme, boolean preprocess, boolean filter, double threshold) throws Exception {
		return metabolizeWithEnzyme(substrate, enzyme, true, preprocess, filter, threshold);		
	}
	
	public ArrayList<Biotransformation> metabolizeWithEnzyme(IAtomContainer substrate, Enzyme enzyme, boolean predictESSpecificity, boolean preprocess, boolean filter, double threshold) throws Exception {
		IAtomContainer clonedSub = substrate.clone();
		
		if(preprocess){
			clonedSub = ChemStructureManipulator.preprocessContainer(clonedSub);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clonedSub);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(clonedSub);	
		} else{
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clonedSub);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(clonedSub);
		}
		if(this.bSystem.getEnzymeHash().containsKey(enzyme.getName())){
			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
			ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(clonedSub);
//			System.err.println("chemClasses :"+ chemClasses);	
			
			if (predictESSpecificity) {
				if(this.esspredictor.isValidSubstrate( clonedSub, enzyme.getName(), chemClasses)){
					biotransformations = this.applyReactionsAndReturnBiotransformations(substrate, enzyme.getReactionSet(), preprocess, filter, threshold);
					for(Biotransformation bt : biotransformations){
						if(bt.getEnzymeNames() == null || bt.getEnzymeNames().isEmpty()){
							ArrayList<String> elist = new ArrayList<String>();
							elist.add(String.valueOf(enzyme.getName()));
							bt.setEnzymeNames(elist);
						} else{
							bt.getEnzymeNames().add(enzyme.getName());
						}
					}				
				} 
				else {
					biotransformations = this.applyReactionsAndReturnBiotransformations(substrate, enzyme.getReactionSet(), preprocess, filter, threshold);
					for(Biotransformation bt : biotransformations){
						if(bt.getEnzymeNames() == null || bt.getEnzymeNames().isEmpty()){
							ArrayList<String> elist = new ArrayList<String>();
							elist.add(String.valueOf(enzyme.getName()));
							bt.setEnzymeNames(elist);
						} 
						else{
							bt.getEnzymeNames().add(enzyme.getName());
						}
					}
				}						
			}			
			return biotransformations;
			
		} else {
			throw new IllegalArgumentException(enzyme.getName() + " is not associated with the biosystem " + this.getBioSystemName());
		}
	}


	public ArrayList<Biotransformation> metabolizeWithEnzyme(IAtomContainer substrate, Enzyme enzyme, boolean predictESSpecificity, boolean preprocess, boolean filter, int nrOfSteps, double threshold) throws Exception {
		IAtomContainerSet substrates = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);		
		substrates.addAtomContainer(substrate);
		IAtomContainerSet containers = substrates;
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		
		int counter = 0;
		
		while(nrOfSteps>0){
			counter++;
			ArrayList<Biotransformation> currentProducts = metabolizeWithEnzyme(containers, enzyme, predictESSpecificity, preprocess, filter, threshold);
			nrOfSteps--;
			if(!currentProducts.isEmpty()){
				biotransformations.addAll(currentProducts);
				containers.removeAllAtomContainers();
				containers = extractProductsFromBiotransformations(currentProducts);
			}
			else {
				break;
			}
		}

		
		return biotransformations;
	}

	public ArrayList<Biotransformation> metabolizeWithEnzyme(IAtomContainerSet substrates, Enzyme enzyme, boolean predictESSpecificity, boolean preprocess, boolean filter, double threshold) throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer atc : substrates.atomContainers()){
			biotransformations.addAll(metabolizeWithEnzyme(atc, enzyme, predictESSpecificity, preprocess, filter, threshold));
		}
		return biotransformations;
	}

	
	public boolean isValidSubstrate(IAtomContainer target, String enzymeName) throws Exception{		
		return esspredictor.isValidSubstrate(target, enzymeName);		
	}
	
	public ArrayList<Biotransformation> metabolizeWithEnzymes(IAtomContainer target,
			ArrayList<Enzyme> enzymes, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		return  metabolizeWithEnzymes(target, enzymes, true, preprocess, filter, scoreThreshold);
	}
	
	public ArrayList<Biotransformation> metabolizeWithEnzymes(IAtomContainer target,
			ArrayList<Enzyme> enzymes, boolean predictESSpecificity, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{			
			ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();
			ArrayList<Enzyme> metabolizingEnzymes = new ArrayList<Enzyme>();
			LinkedHashMap<String, ArrayList<String>> reactToEnzymes = new LinkedHashMap<String, ArrayList<String>>();
			LinkedHashMap<String, MetabolicReaction> reactions = new LinkedHashMap<String, MetabolicReaction>();
			ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
//			IAtomContainer starget = ChemStructureManipulator.standardizeMoleculeWithCopy(target, preprocess);
			IAtomContainer starget = target.clone();
			
			if (preprocess) {
				try {
//					starget = standardizeMoleculeWithCopy(target);
					starget = ChemStructureManipulator.preprocessContainer(target);
					AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
				}
				catch (Exception e) {
					System.out.println(e);
				}
			}
			else{
//				starget = standardizeMoleculeWithCopy(target, false);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(starget);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
				
			}
						
			
			InChIGenerator gen0 = this.inchiGenFactory.getInChIGenerator(target);
			
			if(target.getProperty("InChI") == null || ((String) target.getProperty("InChI")).trim().length()==0){
				target.setProperty("InChI", gen0.getInchi());
				target.setProperty("InChIKey", gen0.getInchiKey());
				Utilities.addPhysicoChemicalProperties(target);
				target.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(target));
			}
			
			if(target.getProperty("SMILES") == null) {
				target.setProperty("SMILES", this.smiGen.create(AtomContainerManipulator.removeHydrogens(target)));
//				System.err.println("target SMILES: " + this.smiGen.create(AtomContainerManipulator.removeHydrogens(target)));			
			}

			
			ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(starget);
			if (predictESSpecificity) {
				for(Enzyme enz : enzymes){
					if(esspredictor.isValidSubstrate(starget, enz.getName(), chemClasses)){
						metabolizingEnzymes.add(enz);
					}	
				}				
			} else {
				metabolizingEnzymes = enzymes;
			}

			
			
			for (Enzyme enzy : metabolizingEnzymes){
				for(MetabolicReaction m : enzy.getReactionSet()){
					if(m.getComonName().equals("EAWAG_RULE_BT1764")) {
					}
					if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, starget)){
						if(reactToEnzymes.get( m.getReactionName() ) == null){
							reactToEnzymes.put(m.getReactionName(), new ArrayList<String>());
							reactToEnzymes.get(m.getReactionName()).add(enzy.getName());
						}
						else{
							reactToEnzymes.get(m.getReactionName()).add(enzy.getName());
						}
						reactions.put(m.getReactionName(), m);	
						matchedReactions.add(m);
					}					
				}			
			}
				
		
			ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();
		
			if(filter == false){
				filteredReactions = matchedReactions;		
			} else{
				filteredReactions = new ArrayList<MetabolicReaction>(this.mRFilter.filterReactions(matchedReactions).values());
			}
			
						
			for(MetabolicReaction j : filteredReactions){
				IAtomContainer n = this.smiParser.parseSmiles(this.smiGen.create(starget));
				IAtomContainerSet partialSet = this.generateAllMetabolitesFromAtomContainer(n, j, true);
				Double score=0.0;
				AtomContainerSet subs = new AtomContainerSet();
				AtomContainerSet prod = new AtomContainerSet();								
				if(partialSet.getAtomContainerCount()>0){
					
					if(target.getProperty("Score") != null){							
						score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(j.name)  );
					}else{
						try{
							score = this.bSystem.getReactionsORatios().get(j.name);
						} catch(Exception e) {
							score = 0.5;
						}
					}
					
					
					if(score>=scoreThreshold){
						subs.addAtomContainer(target);
						for(IAtomContainer pc : partialSet.atomContainers()){
							if(containsK(pc)) continue;
							if(!ChemStructureExplorer.isUnneccessaryMetabolite(pc)){
								try{
									InChIGenerator gen = this.inchiGenFactory.getInChIGenerator(pc);
									pc.setProperty("InChI", gen.getInchi());
									pc.setProperty("InChIKey", gen.getInchiKey());
									pc.setProperty("SMILES", this.smiGen.create(AtomContainerManipulator.removeHydrogens(pc)));
								}catch (CDKException c){
										System.err.println(c.getLocalizedMessage());
								}
								Utilities.addPhysicoChemicalProperties(pc);
								prod.addAtomContainer(AtomContainerManipulator.removeHydrogens(pc));
								prod.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(target));
							}
						}										
						ArrayList<Enzyme> enzList = new ArrayList<Enzyme>();
						Biotransformation bioT = new Biotransformation(subs, j.name, reactToEnzymes.get(j.name), prod, score, this.getBioSystemName());
						results.add(bioT);
					}
				}	
			}
			return results;
	}
	

	public ArrayList<Biotransformation> metabolizeWithEnzymes(IAtomContainerSet substrates, ArrayList<Enzyme> enzymes, boolean preprocess, boolean filter, double threshold) throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer ac : substrates.atomContainers()){
			biotransformations.addAll(metabolizeWithEnzymes(ac, enzymes, preprocess, filter, threshold));
		}
		return biotransformations;
	}	
	
	public ArrayList<Biotransformation> metabolizeWithEnzymesDephtFirst(IAtomContainer target,
			ArrayList<Enzyme> enzymes, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold) throws Exception{
		return metabolizeWithEnzymesDephtFirst(target,
				enzymes, true, preprocess, filter, nrOfSteps, scoreThreshold);
	}
	
	public ArrayList<Biotransformation> metabolizeWithEnzymesDephtFirst(IAtomContainer target,
		ArrayList<Enzyme> enzymes, boolean predictESSpecificity, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();
		for(Enzyme enz : enzymes){
			results.addAll(metabolizeWithEnzyme(target, enz, predictESSpecificity, preprocess, filter, nrOfSteps, scoreThreshold) );		
		}
		return results;
	}	
	
	public ArrayList<Biotransformation> metabolizeWithEnzymesBreadthFirst(IAtomContainer target, ArrayList<Enzyme> enzymes, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold) throws Exception{
		IAtomContainerSet targets = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		targets.addAtomContainer(target);
		return  metabolizeWithEnzymesBreadthFirst(targets, enzymes, preprocess, filter, nrOfSteps, scoreThreshold);		
	}
	
	public ArrayList<Biotransformation> metabolizeWithEnzymesBreadthFirst(IAtomContainerSet targets,
			ArrayList<Enzyme> enzymes, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold) throws Exception{
			ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();
			IAtomContainerSet containers = targets;

			int counter = 0;
			
			while(nrOfSteps>0){
				counter++;
				ArrayList<Biotransformation> currentProducts = metabolizeWithEnzymes(containers, enzymes, preprocess, filter, scoreThreshold);
				nrOfSteps--;
				
				if(!currentProducts.isEmpty()){
					results.addAll(currentProducts);
					containers.removeAllAtomContainers();
					containers = extractProductsFromBiotransformations(currentProducts);
				}
				else{
					break;
				}				
			}
			return results;
		}
	
	
	public ArrayList<Biotransformation> applyReactionAtOnceAndReturnBiotransformations(IAtomContainer target,
			MetabolicReaction reaction, boolean preprocess, Double scoreThreshold) throws Exception{
		
//		System.err.println(reaction.getReactionName());
		ArrayList<Biotransformation> results = new ArrayList<Biotransformation>();		
		IAtomContainer starget = target.clone();
//		IAtomContainer starget = this.standardizeMoleculeWithCopy(target);
		
		if (preprocess) {
			try {
				starget = ChemStructureManipulator.preprocessContainer(starget);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
			}
			catch (Exception e) {
				System.out.println(e);
			}
		}
		else{
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(starget);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(starget);
		}
//		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
//		ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();		

		InChIGenerator gen0 = this.inchiGenFactory.getInChIGenerator(target);
		target.setProperty("InChI", gen0.getInchi());
		target.setProperty("InChIKey", gen0.getInchiKey());
		target.setProperty("SMILES", this.smiGen.create(target));
		target.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(target));
				
		IAtomContainerSet partialSet = generateAllMetabolitesFromAtomContainerViaTransformationAtAllLocations(
				starget, reaction.getSmirksReaction(), false);
		
		
		Double score=0.0;
		IAtomContainerSet subs = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet prods = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		
		if(partialSet.getAtomContainerCount()>0){
			
			if(target.getProperty("Score") !=null){	
				
				score = ( new Double((Double) target.getProperty("Score")) * this.bSystem.getReactionsORatios().get(reaction.name)  );
			}else{
				score = this.bSystem.getReactionsORatios().get(reaction.name);
			}
			if(score>=scoreThreshold){
				subs.addAtomContainer(target);
				for(IAtomContainer pc : partialSet.atomContainers()){
					
					InChIGenerator gen = this.inchiGenFactory.getInChIGenerator(pc);
					pc.setProperty("InChI", gen.getInchi());
					pc.setProperty("InChIKey", gen.getInchiKey());
					pc.setProperty("SMILES", this.smiGen.create(pc));
					pc.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(pc));
					prods.addAtomContainer(AtomContainerManipulator.removeHydrogens(pc));
				}
				
				Biotransformation bioT = new Biotransformation(subs, reaction.name, null, prods, score, this.getBioSystemName() );
				results.add(bioT);
			}
		}
		return results;
	}


	public ArrayList<Biotransformation> applyReactionAtOnceAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();	
		
		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
		ArrayList<MetabolicReaction> filteredReactions = new ArrayList<MetabolicReaction>();
		
		for (MetabolicReaction i : reactions) {
			boolean match_constraints = ChemStructureExplorer.compoundMatchesReactionConstraints(i, target);
			if (match_constraints) {
				
				matchedReactions.add(i);
			}
		}		
		
		if(filter == false){
			filteredReactions = matchedReactions;		
		} 
		else{
			filteredReactions = new ArrayList<MetabolicReaction>(this.mRFilter.filterReactions(matchedReactions).values());
		}
		
		for(MetabolicReaction m : matchedReactions) {
			//System.out.println(m.getReactionName());
		}
		//System.out.println("matchedReactions : " + matchedReactions.size());
		
		for(MetabolicReaction f : filteredReactions) {
			//System.out.println(f.getReactionName());
		}
		//System.out.println("filteredReactions : " + filteredReactions.size());

		for(MetabolicReaction mreact : reactions) {
			
			ArrayList<Biotransformation> bt = applyReactionAtOnceAndReturnBiotransformations(target, mreact, preprocess, scoreThreshold);
//			System.err.println("bt " + bt.get(0).getProducts().getAtomContainerCount());
			ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();
			
			for(Biotransformation b : bt){
				if(b.getSubstrates().getAtomContainerCount() == b.getProducts().getAtomContainerCount() &&
						b.getSubstrates().getAtomContainerCount() == 1 &&
						ChemStructureExplorer.inchikeyEqualityHolds(b.getSubstrates().getAtomContainer(0), b.getProducts().getAtomContainer(0))){
					
					System.err.println("Removing " + b.getReactionType());
				} else{
					selectedBiotransformations.add(b);
				}
			}	
			
			biotransformations.addAll(selectedBiotransformations);			
		}
			
		
		return biotransformations;	
	}
	
	public ArrayList<Biotransformation> applyReactionAtOnceAndReturnBiotransformations(IAtomContainerSet targets,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
        	
		for(IAtomContainer a : targets.atomContainers()) {
			biotransformations.addAll(applyReactionAtOnceAndReturnBiotransformations(a, reactions, preprocess, filter, scoreThreshold));			
		}
		return biotransformations;		
	}
	
	public ArrayList<Biotransformation> applyReactionAtOnceAndReturnBiotransformations(IAtomContainer target,
			ArrayList<MetabolicReaction> reactions, boolean preprocess, boolean filter, int nrOfSteps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		int step = 0;
		int i = nrOfSteps;		 
		IAtomContainerSet startingSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		startingSet.addAtomContainer(target);		
		while(nrOfSteps>0) {			
			step++;
			System.out.println("Step " + step + " out of " + i);
			ArrayList<Biotransformation> partialBiotransf = applyReactionAtOnceAndReturnBiotransformations(target,
					reactions, preprocess, filter, scoreThreshold);
			nrOfSteps--;
			System.out.println("Remaining steps " + nrOfSteps);
			if(partialBiotransf.size()>0){
				biotransformations.addAll(partialBiotransf);
				startingSet.removeAllAtomContainers();
				startingSet = this.extractProductsFromBiotransformations(partialBiotransf);
				
				for(IAtomContainer a : startingSet.atomContainers()){
//					biotransformations.addAll(applyReactionAtOnceAndReturnBiotransformations(standardizeMoleculeWithCopy(a),
//							reactions, false, filter, nrOfSteps, scoreThreshold));
					biotransformations.addAll(applyReactionAtOnceAndReturnBiotransformations(a,
							reactions, false, filter, nrOfSteps, scoreThreshold));
				}

				
			} else {
				break;
			}
			
		}
		
		return biotransformations;
		
	}

	public SmilesParser getSmiParser() {
		return smiParser;
	}

	public IAtomContainerSet extractProductsFromBiotransformations(ArrayList<Biotransformation> biotransformations) throws Exception{
		IAtomContainerSet acontainers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		LinkedHashMap<String, IAtomContainer> hMap = new LinkedHashMap<String, IAtomContainer>();
		for(Biotransformation b : biotransformations){
			int numAtoms_substrate = 0;
			IAtomContainerSet substrates = b.getSubstrates();
			for(int i = 0; i < substrates.getAtomContainerCount(); i++){
				IAtomContainer oneSubstrate = substrates.getAtomContainer(i);
				//AtomContainerManipulator.suppressHydrogens(oneSubstrate);
				//numAtoms_substrate += oneSubstrate.getAtomCount(); 
				numAtoms_substrate += countNonHydrogenAtoms(oneSubstrate);
			}
			for(IAtomContainer ac : b.getProducts().atomContainers()){
				String ikey = ac.getProperty("InChIKey");			
				//The suppressHydrogens function may not work for ac. The assumption here is not that useful and doesn't make sense
				//AtomContainerManipulator.suppressHydrogens(ac);
				//int numAtoms_metabolites = ac.getAtomCount();
				int numAtoms_metabolites = countNonHydrogenAtoms(ac);
				//Use 5 here considering the number of non-hydrogen atoms in sulfate (was glycine)
				if(numAtoms_metabolites >= numAtoms_substrate + 4 || ChemStructureExplorer.getMajorIsotopeMass(ac) > 1500){
					ac.setProperty("isEndProduct", true);
				}
				if(!hMap.containsKey(ikey)){
					hMap.put(ikey, ac);
					
				}
			}
		}
		ArrayList<IAtomContainer> at =  new ArrayList<IAtomContainer>( hMap.values());
//		System.err.println(at.get(0).getClass());
		for(IAtomContainer a : at){
			if(countNonHydrogenAtoms(a) <= 4) continue;
			acontainers.addAtomContainer(ChemStructureManipulator.preprocessContainer(a));
		}
		return acontainers;
//		return ChemStructureExplorer.uniquefy(acontainers);
	}
	/**
	 * This 
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public int countNonHydrogenAtoms(IAtomContainer molecule) throws Exception{
		int count = 0;
		for(int i = 0; i < molecule.getAtomCount(); i++){
			if(!molecule.getAtom(i).getSymbol().equalsIgnoreCase("H")) count++;
		}
		return count;
	}
	
	public IAtomContainerSet extractProductsFromBiotransformationsWithTransformationData(ArrayList<Biotransformation> biotransformations, LinkedHashMap<String, MetabolicReaction> customReactionHash, boolean annotate) throws Exception{
		ArrayList<Biotransformation> uniqueBiotransformations = Utilities.selectUniqueBiotransformations(biotransformations);
		IAtomContainerSet acontainers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		LinkedHashMap<String, IAtomContainer> hMap = new LinkedHashMap<String, IAtomContainer>();
		//System.out.println("Biotransformations: " + biotransformations.size());
		System.out.println("Unique Biotransformations: " + uniqueBiotransformations.size());
		System.out.println("Unique metabolites: " + Utilities.extractProductsFromBiotransformations(uniqueBiotransformations).getAtomContainerCount());
		int metaboliteID = 0;
		if(uniqueBiotransformations != null){
			for(Biotransformation b : uniqueBiotransformations){
				for(IAtomContainer ac : b.getProducts().atomContainers()){
					if(containsK(ac)) continue;
					IAtomContainer hash_ac;
					LinkedHashMap<Object, Object> properties = new LinkedHashMap<Object, Object>();
					String ikey = ac.getProperty("InChIKey");
					
					if(hMap.containsKey(ikey)) {
						hash_ac = hMap.get(ikey).clone();
						properties.put("InChI", hash_ac.getProperty("InChI"));
						properties.put("InChIKey", hash_ac.getProperty("InChIKey"));						
						properties.put("SMILES", hash_ac.getProperty("SMILES"));							
						properties.put("Synonyms", hash_ac.getProperty("Synonyms"));
						properties.put("PUBCHEM_CID", hash_ac.getProperty("PUBCHEM_CID"));
						properties.put("Molecular formula", hash_ac.getProperty("Molecular formula"));	
						properties.put("Major Isotope Mass", hash_ac.getProperty("Major Isotope Mass"));
						properties.put("ALogP", hash_ac.getProperty("ALogP"));
						properties.put("Lipinski_Violations", hash_ac.getProperty("Lipinski_Violations"));
						properties.put("Insecticide_Likeness_Violations", hash_ac.getProperty("Insecticide_Likeness_Violations"));
						properties.put("Post_Em_Herbicide_Likeness_Violations", hash_ac.getProperty("Post_Em_Herbicide_Likeness_Violations"));						
						properties.put("Metabolite ID", hash_ac.getProperty("Metabolite ID"));
						properties.put(CDKConstants.TITLE, hash_ac.getProperty(CDKConstants.TITLE));						
					}
					else {
						hash_ac = ac.clone();
						LinkedHashMap<Object, Object> refProps = new LinkedHashMap<Object, Object>();
						String synonyms = null;
						String pubchemCID = null;
								
						if(annotate){
							
							System.out.println("\n\n===========================================");
							System.out.println("Fetching CIDs and synonyms from PubChem");
							System.out.println("===========================================\n\n");
							
							LinkedHashMap<String,ArrayList<String>> data = ChemdbRest.getSynonymsObjectViaInChIKey(ikey);
							
							if(data != null ){
								pubchemCID 	=  data.get("CID").get(0);
								synonyms 	= StringUtils.join(data.get("Synonyms"), "\n");
							}
						}						
						
						refProps.put("InChI", ac.getProperty("InChI"));
						refProps.put("InChIKey", ac.getProperty("InChIKey"));
						refProps.put("SMILES", hash_ac.getProperty("SMILES"));
						refProps.put("Synonyms", synonyms);
						refProps.put("PUBCHEM_CID", pubchemCID);
						
						if(ac.getProperty("Molecular formula") !=null){
							refProps.put("Molecular formula", ac.getProperty("Molecular formula"));	
						}
						else
						{
							refProps.put("Molecular formula", ChemStructureExplorer.getMolecularFormula(ac));	
						}

						if(ac.getProperty("Molecular formula") !=null){
							refProps.put("Major Isotope Mass", ac.getProperty("Major Isotope Mass"));						
							refProps.put("ALogP", ac.getProperty("ALogP"));									
							refProps.put("Lipinski_Violations", ac.getProperty("Lipinski_Violations"));
							refProps.put("Insecticide_Likeness_Violations", ac.getProperty("Insecticide_Likeness_Violations"));
							refProps.put("Post_Em_Herbicide_Likeness_Violations", ac.getProperty("Post_Em_Herbicide_Likeness_Violations"));

						}
						else {
							LinkedHashMap<String, String> physchem = ChemStructureExplorer.computePhysicoChemicalProperties(ac);
							refProps.put("Major Isotope Mass", physchem.get("Major Isotope Mass"));	
							refProps.put("ALogP", physchem.get("ALogP"));	

							LinkedHashMap<String, Integer> violations = ChemStructureExplorer.calculateLikenessViolations(ac);
							refProps.put("Lipinski_Violations", String.format("%s", violations.get("Lipinski_Violations")));
							refProps.put("Insecticide_Likeness_Violations", String.format("%s", violations.get("Insecticide_Likeness_Violations")));
							refProps.put("Post_Em_Herbicide_Likeness_Violations", String.format("%s", violations.get("Post_Em_Herbicide_Likeness_Violations")));

						}
							
						metaboliteID++;
						refProps.put("Metabolite ID", "BTM" + String.format("%05d", metaboliteID));
						refProps.put(CDKConstants.TITLE, "BTM" + String.format("%05d", metaboliteID));
						hash_ac.setProperties(refProps);
						hMap.put((String) ac.getProperty("InChIKey"), hash_ac);
						
						properties = (LinkedHashMap<Object, Object>) refProps.clone();				
					}
					/**
					 * Edited by Siyang
					 */
					//String reactionName = customReactionHash.get(b.getReactionType().toString()).getComonName();
					String reactionName;
					if(customReactionHash.get(b.getReactionType().toString())!=null){
						reactionName = customReactionHash.get(b.getReactionType().toString()).getComonName();
					}
					else {
						if(b.getReactionType().toString().contains(" AndFromCyProduct")) {
							String reactionType = b.getReactionType().toString().replace(" AndFromCyProduct", "");
							if(customReactionHash.get(reactionType)!=null) {
								reactionName = customReactionHash.get(reactionType).getComonName();
							}
							else reactionName = b.getReactionType();
							reactionName  = reactionName + " AndFromCyProduct";
						}
						else reactionName = b.getReactionType();
					}
					
					if (reactionName.length() == 0){
						reactionName = customReactionHash.get(b.getReactionType().toString()).toString();
					}
					properties.put("Reaction",reactionName);
					
					if(customReactionHash.get(b.getReactionType().toString())!=null){
						properties.put("Reaction ID",customReactionHash.get(b.getReactionType().toString()).getBTRMID());
					}
					else{
						if(ac.getProperty("Reaction ID")!=null) properties.put("Reaction ID",ac.getProperty("Reaction ID"));
						else properties.put("Reaction ID","N/A");
					}					
					if(b.getEnzymeNames().size()>0){
						ArrayList<String> enzymes = new ArrayList<String>();										
						for(String en : b.getEnzymeNames()){
							if(en.toString().contains("EC_")){
								enzymes.add(en.toString().replace("EC_", "EC ").replaceAll("_", ".")) ;
							}
							else {
		 						String enzymeName = en.toString().
									replace("HYDROXYCINNAMATE_DECARBOXYLASE", "4-Hydroxycinnamate decarboxylase").								
									replace("PHENOLIC_ACID_DECARBOXYLASE", "Bacterial phenolic acid decarboxylase (EC 4.1.1.-)").
									replace("DECARBOXYLASE", "Unspecified bacterial decarboxylase").
									replace("DEMETHYLASE", "Unspecified bacterial demethylase").
									replace("DEHYDROXYLASE", "Unspecified bacterial dehydroxylase").
									replace("DEHYDROXYLASE", "Unspecified bacterial dehydroxylase").
									replace("BACTERIAL_LACTONASE", "Unspecified bacterial lactonase").					
									replace("VINYL_PHENOL_REDUCTASE", "Vinyl phenol reductase").
									replace("ABKAR1", "Alpha,beta-ketoalkene double bond reductase 1").
									replace("UDP_GLUCURONOSYLTRANSFERASE", "Bacterial UDP-glucuronosyltransferase").
									replace("ACETYLTRANSFERASE", "Unspecified acetyltransferase").
									replace("UNSPECIFIED_BACTERIAL_ISOFLAVONE_REDUCTASE", "Unspecified bacterial isoflavone reductase").
									replace("UNSPECIFIED_GUT_BACTERIAL_ENZYME", "Unspecified gut bacterial enzyme").
									replace("UNSPECIFIED_ENVIRONMENTAL_BACTERIAL_ENZYME", "Unspecified environmental bacterial enzyme").
									replace("BACTERIAL_BILE_SALT_3_HYDROXYSTEROID_DEHYDROGENASE", "Bacterial bile salt 3-hydroxysteroid dehydrogenase").
									replace("BACTERIAL_BILE_SALT_7_HYDROXYSTEROID_DEHYDROGENASE", "Bacterial bile salt 7-hydroxysteroid dehydrogenase").
									replace("BACTERIAL_BILE_SALT_12_HYDROXYSTEROID_DEHYDROGENASE", "Bacterial bile salt 12-hydroxysteroid dehydrogenase").
									replace("BACTERIAL_NITROREDUCTASE", "Bacterial oxygen-insensitive NADPH nitroreductase").
									replace("EC_3_5_1_24", "Choloylglycine hydrolase");
										
								enzymes.add(enzymeName.toString());
							}
						}
		//				ac.setProperty("Enzyme(s)", StringUtils.join(enzymes, "\n"));
						properties.put("Enzyme(s)", StringUtils.join(enzymes, "\n"));
						properties.put("Biosystem", b.getBioSystemName().name());
					}
					if(b.getSubstrates().getAtomContainerCount() == 1){
		//				System.err.println(b.getSubstrates().getAtomContainer(0));
						
						IAtomContainer substrate = hMap.get(b.getSubstrates().getAtomContainer(0).getProperty("InChIKey"));
						String tt = null;
						
						if (substrate == null){
							substrate = b.getSubstrates().getAtomContainer(0).clone();
							LinkedHashMap<Object, Object> refProps = new LinkedHashMap<Object, Object>();
							
							if(substrate.getProperty("SMILES") != null) {
								refProps.put("SMILES", substrate.getProperty("SMILES"));
							}
							else {								
								refProps.put("SMILES", this.smiGen.create(substrate));
							}							
							
							
							if(substrate.getProperty("Molecular formula") !=null){
								refProps.put("Molecular formula", ac.getProperty("Molecular formula"));	
							}
							else
							{
								refProps.put("Molecular formula", ChemStructureExplorer.getMolecularFormula(substrate));	
							}
							
							if(substrate.getProperty("Major Isotope Mass") !=null){
								refProps.put("Major Isotope Mass", substrate.getProperty("Major Isotope Mass"));						
								refProps.put("ALogP", substrate.getProperty("ALogP"));
								refProps.put("Lipinski_Violations", substrate.getProperty("Lipinski_Violations"));
								refProps.put("Insecticide_Likeness_Violations", substrate.getProperty("Insecticide_Likeness_Violations"));
								refProps.put("Post_Em_Herbicide_Likeness_Violations", substrate.getProperty("Post_Em_Herbicide_Likeness_Violations"));

							}
							else {
								LinkedHashMap<String, String> physchem = ChemStructureExplorer.computePhysicoChemicalProperties(substrate);							
								refProps.put("Major Isotope Mass", String.format("%.8s", Double.valueOf(physchem.get("Major Isotope Mass"))));	
								refProps.put("Major Isotope Mass", String.format("%.8s", Double.valueOf(physchem.get("Major Isotope Mass"))));	
								
								LinkedHashMap<String, Integer> violations = ChemStructureExplorer.calculateLikenessViolations(substrate);
								refProps.put("Lipinski_Violations", String.format("%s", violations.get("Lipinski_Violations")));
								refProps.put("Insecticide_Likeness_Violations", String.format("%s", violations.get("Insecticide_Likeness_Violations")));
								refProps.put("Post_Em_Herbicide_Likeness_Violations", String.format("%s", violations.get("Post_Em_Herbicide_Likeness_Violations")));

							
							}					
						
							tt = (String) substrate.getProperty(CDKConstants.TITLE);
							if (tt == null){
								String syno = (String) substrate.getProperty("Synonyms");
								if(syno != null){
									tt = Utilities.returnFirstCleanSynonym(syno.split("\n"));											
								}
								if(tt == null){
									tt = (String) substrate.getProperty("Name");				
									if(tt == null) {
										
										tt = (String) substrate.getProperty("NAME");
										
										if (tt == null ){
											tt = (String) substrate.getProperty("DATABASE_ID");
											
											if(tt == null) {
												tt = (String) substrate.getProperty("DRUGBANK_ID");
			
												if(tt == null){
													tt = (String) substrate.getProperty("Metabolite ID");
			
													if(tt == null){
														tt = (String) substrate.getProperty("$MolName");
														
														if (annotate && tt == null){
															System.out.println("\n\n===========================================");
															System.out.println("Fetching CIDs and synonyms from PubChem");
															System.out.println("===========================================\n\n");
															
															LinkedHashMap<String,ArrayList<String>> data = ChemdbRest.getSynonymsObjectViaInChIKey(ikey);
															
															if(data != null ){
//																String pubchemCID_ 	=  data.get("CID").get(0);
//																String synonyms		= StringUtils.join(data.get("Synonyms"), "\n");
																if (data.get("Synonyms") != null){
																	tt = Utilities.returnFirstCleanSynonym(data.get("Synonyms"));
//																	System.err.println("SYNO TT 2 IS: " + tt);
																}
															}
														}
													}
												}
											}
																			
										}							
									}
									
								}			
							}
							refProps.put(CDKConstants.TITLE, tt);
							substrate.addProperties(refProps);
							hMap.put((String) substrate.getProperty("InChIKey"), substrate);
						}
						properties.put("Precursor ID", substrate.getProperty(CDKConstants.TITLE));	
						properties.put("Precursor SMILES", substrate.getProperty("SMILES"));
						properties.put("Precursor InChI", substrate.getProperty("InChI"));
						properties.put("Precursor InChIKey", substrate.getProperty("InChIKey"));
						properties.put("Precursor ALogP", substrate.getProperty("ALogP"));
						properties.put("Precursor Major Isotope Mass", substrate.getProperty("Major Isotope Mass"));				
					}
	//				
					hash_ac.setProperties(properties);
					acontainers.addAtomContainer(hash_ac);
				
				}
			}
		}	
		return acontainers;
		
	}

	public IAtomContainerSet extractProductsFromBiotransformationsWithTransformationData(ArrayList<Biotransformation> biotransformations, boolean annotate) throws Exception{
		return  extractProductsFromBiotransformationsWithTransformationData(biotransformations, this.reactionsHash, annotate);
	}
	
	

	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName) throws Exception{
		saveBioTransformationProductsToSdf(biotransformations, outputFileName, this.reactionsHash, false);

	}
	
	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception{
		saveBioTransformationProductsToSdf(biotransformations, outputFileName, this.reactionsHash, annotate);

	}

	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName,  LinkedHashMap<String, MetabolicReaction> customReactionHash) throws Exception{
		saveBioTransformationProductsToSdf(biotransformations, outputFileName, customReactionHash, false);

	}
	
	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName, LinkedHashMap<String, MetabolicReaction> customReactionHash, boolean annotate) throws Exception{
		IAtomContainerSet uniqueSetOfProducts = extractProductsFromBiotransformationsWithTransformationData(biotransformations, customReactionHash, annotate);
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));		
		sdfWriter.write(uniqueSetOfProducts);
		sdfWriter.close();
	}
	
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformationChain(IAtomContainer target,
			String pathwayName, boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		
		return applyPathwaySpecificBiotransformationsChain(target, pathwayName, preprocess, filter,  nr_of_steps, 0.0);

	}
	
	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName) throws Exception{
		saveBioTransformationProductsToCSV(biotransformations, outputFileName, this.reactionsHash, false);

	}
	
	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception{
		saveBioTransformationProductsToCSV(biotransformations, outputFileName, this.reactionsHash, annotate);

	}

	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName,  LinkedHashMap<String, MetabolicReaction> customReactionHash) throws Exception{
		saveBioTransformationProductsToSdf(biotransformations, outputFileName, customReactionHash, false);

	}
	
	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName, LinkedHashMap<String, MetabolicReaction> customReactionHash, boolean annotate) throws Exception{
		try{
			IAtomContainerSet products = extractProductsFromBiotransformationsWithTransformationData(biotransformations, customReactionHash, annotate);
			FileUtilities.saveAtomContainerSetToCSV(products, outputFileName);			
		}
		catch (Exception e) {
			e.printStackTrace();
		}

			
	}
	
	public void saveBioTransformationProductsToSDF(ArrayList<Biotransformation> biotransformations, String outputFileName, LinkedHashMap<String, MetabolicReaction> customReactionHash, boolean annotate) throws Exception{
		try{
			IAtomContainerSet products = extractProductsFromBiotransformationsWithTransformationData(biotransformations, customReactionHash, annotate);
			FileUtilities.saveAtomContainerSetToSDF(products, outputFileName);	
		}
		catch (Exception e) {
			e.printStackTrace();
		}

			
	}
	
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformations(IAtomContainer target,
	String pathwayName, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		return applyPathwaySpecificBiotransformations(target,
				pathwayName, true, preprocess, filter, scoreThreshold);
	}
	
	/**
	 * Predicts the metabolism of a compound for a specific pathway, and returns metabolites with a minimal given score threshold.
	 * @param target - An AtomContainer that represent a chemical compound
	 * @param pathwayName - the name of a MetabolicPathway object.
	 * @param preprocess - specifies whether the compounds must be pre-processed.
	 * @param filter - specifies whether the reactions should be filtered, according to priority rules.
	 * @param scoreThreshold - the minimal score that a metabolite must have to be returned.
	 * @return an ArrayList of biotransformations
	 * @throws Exception
	 *  			  : Throws an Exception
	 */
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformations(IAtomContainer target, String pathwayName, boolean predictESSpecificity, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();				
		if(ChemStructureExplorer.isBioTransformerValid(target)){		
			ArrayList<Enzyme> enzymes = this.bSystem.getMetPathwaysHash().get(pathwayName);
			for(Enzyme e : enzymes){
				biotransformations.addAll(this.metabolizeWithEnzyme(target, e, predictESSpecificity, preprocess, filter, scoreThreshold));
			}			
		}

		return biotransformations;
	}
	
	/**
	 * Predicts the metabolism of a compound for a specific pathway, and returns metabolites with a minimal given score threshold.
	 * @param targets - A set of chemical compounds
	 * @param pathwayName - the name of a MetabolicPathway object.
	 * @param preprocess - specifies whether the compounds must be pre-processed.
	 * @param filter - specifies whether the reactions should be filtered, according to priority rules.
	 * @param scoreThreshold - the minimal score that a metabolite must have to be returned.
	 * @return an ArrayList of biotransformations
	 * @throws Exception
	 *  			  : Throws an Exception
	 */	
	
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformations(IAtomContainerSet targets,
	String pathwayName, boolean preprocess, boolean filter, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet readyTargets = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		for(IAtomContainer atc : targets.atomContainers()){
			IAtomContainer a = atc.clone();
			if(preprocess){
				a = ChemStructureManipulator.preprocessContainer(a);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(a);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(a);	
				readyTargets.addAtomContainer(a);
			} else{
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(a);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(a);
				readyTargets.addAtomContainer(a);
			}
		}
		
		
		for(IAtomContainer aa : readyTargets.atomContainers()){
//			System.err.println("preprocessed target again: " + this.smiGen.create(aa));
			biotransformations.addAll(this.applyPathwaySpecificBiotransformations(aa, pathwayName, false, filter, scoreThreshold));
		}
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformationsChain(IAtomContainerSet targets,
			String pathwayName, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = targets;
		for(IAtomContainer atc : targets.atomContainers()){
			int counter = 0;
			while(nr_of_steps>counter){
				counter++;			
				ArrayList<Biotransformation> currentProducts = applyPathwaySpecificBiotransformations(containers, pathwayName, preprocess, filter, scoreThreshold);
//				nr_of_steps--;
//				System.err.println(currentProducts.size() + " biotransformations at step " + counter);
				if(!currentProducts.isEmpty()){		
					biotransformations.addAll(currentProducts);
					containers.removeAllAtomContainers();
					containers = extractProductsFromBiotransformations(currentProducts);
				}
				else{
					break;
				}
			}
			return biotransformations;			
		}
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
			
	
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformationsChain(IAtomContainer target,
			String pathwayName, boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet targets = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		targets.addAtomContainer(target);
		biotransformations = applyPathwaySpecificBiotransformationsChain(targets,
				pathwayName, preprocess, filter,nr_of_steps, scoreThreshold);
		return Utilities.selectUniqueBiotransformations(biotransformations);
		
	}
	
	
	/**
	 * Here we don't want to consider the metabolites that have too few atoms in the structure.
	 * @param oneMole
	 * @return
	 */
	public boolean isValidMetabolte(IAtomContainer oneMole) throws Exception{
		
		if(this.invalidSMARTS.containInvalidSubstructure(oneMole)) return false;
		int countAtom_nonH = 0;
		int countCarbon = 0;
		for(int t = 0; t < oneMole.getAtomCount(); t++){
			if(!oneMole.getAtom(t).getSymbol().equalsIgnoreCase("H")) countAtom_nonH++;
			if(oneMole.getAtom(t).getSymbol().equalsIgnoreCase("C")) countCarbon++;
		}
		if(countAtom_nonH >=1 && countCarbon >=1) return true;
		else return false;
	}
	/**
	 * This function will always check if each molecule in the molecules set is a polymer, and if so, convert each polymer to monomers.
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> processPolymer(IAtomContainerSet molecules) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<>();
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer startingCompound = molecules.getAtomContainer(i).clone();
			if(this.hp.isPolymer(startingCompound)){
				//substratePool.add(this.hp.convertPolymerToMonomer(startingCompound));
				ArrayList<Biotransformation> monomers_biotransformation = this.hp.converPolymerToMonomer_biotransformation(startingCompound);
				biotransformations.addAll(monomers_biotransformation);
			}
		}
		return biotransformations;
	}
	
	public ArrayList<Biotransformation> processPolymer_oneMole(IAtomContainer molecule) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<>();
		if(this.hp.isPolymer(molecule)){
			//substratePool.add(this.hp.convertPolymerToMonomer(startingCompound));
			ArrayList<Biotransformation> monomers_biotransformation = this.hp.converPolymerToMonomer_biotransformation(molecule);
			biotransformations.addAll(monomers_biotransformation);
		}
		return biotransformations;
	}
	
	public boolean containsK(IAtomContainer ac)  throws Exception{
		for(int t = 0; t < ac.getAtomCount(); t++) {
			if(ac.getAtom(t).getSymbol().equalsIgnoreCase("K")) {
				return true;
			}
		}
		return false;
	}
	/**
	 * Added by Siyang
	 * Convert the cyProduct predicted Results to BioTransformations
	 * @param containers
	 * @param nrOfSteps
	 * @param scoreThreshold
	 * @param outputFolder
	 * @param annotate
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> convertMLProductToBioTransformation(IAtomContainerSet substrates, IAtomContainerSet cyProductMolecules) throws Exception{
		//Biotransformation(IAtomContainerSet substrates, String reactionType, ArrayList<String> enzymeNames, 
		//IAtomContainerSet products, BioSystemName bsysName)
		IAtomContainerSet cleanedCyProductMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class); 
		ArrayList<String> checkedInChIKeyList = new ArrayList<>();
		for(int i = 0; i < cyProductMolecules.getAtomContainerCount(); i++){
			IAtomContainer oneMetabolite = cyProductMolecules.getAtomContainer(i);
			if(!checkedInChIKeyList.contains(oneMetabolite.getProperty("InChIKey"))) {
				checkedInChIKeyList.add(oneMetabolite.getProperty("InChIKey"));
				cleanedCyProductMolecules.addAtomContainer(oneMetabolite);
			}
		}
		ArrayList<Biotransformation> convertedBioTrans = new ArrayList<>();
		HashMap<String, IAtomContainerSet> reactionEnzyme_Metabolite_map = new HashMap<>();
		for(int i = 0; i < cleanedCyProductMolecules.getAtomContainerCount(); i++){
			IAtomContainer oneMolecule = cyProductMolecules.getAtomContainer(i);
			String reactionType = oneMolecule.getProperty("ReactionType");
			String enzymeNameString = oneMolecule.getProperty("Enzyme");
			String key = reactionType + ";" + enzymeNameString;
			if(!reactionEnzyme_Metabolite_map.containsKey(key)) {
				IAtomContainerSet notIncludedMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
				notIncludedMolecules.addAtomContainer(oneMolecule);
				reactionEnzyme_Metabolite_map.put(key, notIncludedMolecules);
			}
			else {
				IAtomContainerSet notIncludedMolecules = reactionEnzyme_Metabolite_map.get(key);
				notIncludedMolecules.addAtomContainer(oneMolecule);
				reactionEnzyme_Metabolite_map.put(key, notIncludedMolecules);
			}
		}
		for(String reactionEnzyme : reactionEnzyme_Metabolite_map.keySet()) {
			IAtomContainerSet notIncludedMolecules = reactionEnzyme_Metabolite_map.get(reactionEnzyme);
			String reactionType = reactionEnzyme.split(";")[0];
			String enzymeNameString = reactionEnzyme.split(";")[1];
			String[] enzymeArray = enzymeNameString.split(" ");
			ArrayList<String> enzymeNames = new ArrayList<>();
			for(int j = 0; j < enzymeArray.length; j++){
				enzymeNames.add(enzymeArray[j]);
			}
			BioSystemName bioSys = BioSystemName.HUMAN;
			Biotransformation bioTrans = new Biotransformation(substrates, reactionType, enzymeNames, notIncludedMolecules, bioSys);
			convertedBioTrans.add(bioTrans);
		}
		return convertedBioTrans;
	}
	/**
	 * Added by Siyang
	 * Find the metabolites in set A and not in set B
	 * The input are two IAtomcontainerSets, one is set A and the other is set B
	 * @param containers
	 * @param nrOfSteps
	 * @param scoreThreshold
	 * @param outputFolder
	 * @param annotate
	 * @throws Exception
	 */
	public IAtomContainerSet notIncludedMolecules(IAtomContainerSet mlPredictedResults, IAtomContainerSet ruleBasedResults) throws Exception{
		IAtomContainerSet notIncludedMolecuels = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i <mlPredictedResults.getAtomContainerCount(); i++){
			if(!ChemStructureExplorer.atomContainerInclusionHolds(ruleBasedResults, mlPredictedResults.getAtomContainer(i))){
				notIncludedMolecuels.addAtomContainer(mlPredictedResults.getAtomContainer(i));
			}
		}
		return notIncludedMolecuels;
		
	}
}
