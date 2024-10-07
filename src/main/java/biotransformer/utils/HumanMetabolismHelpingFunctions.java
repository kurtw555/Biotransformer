package biotransformer.utils;

import java.util.ArrayList;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.dbrelevant.RetriveFromDB;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemicalClassFinder.ChemicalClassName;
import biotransformer.validateModels.InValidSMARTS;
import biotransformer.validateModels.IsEndProductSMARTS;
import biotransformer.validateModels.MultivariateDistribution;

public class HumanMetabolismHelpingFunctions {
	boolean useDB;
	public SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
	public SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	public Cyp450BTransformer cyb;// = new Cyp450BTransformer(BioSystemName.HUMAN);
	public Phase2BTransformer p2b;// = new Phase2BTransformer(BioSystemName.HUMAN);
	public ECBasedBTransformer ecb;//	= new ECBasedBTransformer(BioSystemName.HUMAN);
	public HGutBTransformer hgb;// = new HGutBTransformer();
	HandlePolymers hp = new HandlePolymers();
	RetriveFromDB rfdb;
	InValidSMARTS invalidSMARTS = new InValidSMARTS();
	IsEndProductSMARTS isEndProductSMARTS = new IsEndProductSMARTS();
	//MultivariateDistribution mnd_Liver;
	//MultivariateDistribution mnd_ECBased;
	//MultivariateDistribution mnd_GutSubstrate;
	//MultivariateDistribution mnd_GutMetabolite;
	int cyp450Mode = 1;
	Double scoreThreshold = 0.0;
	boolean useSubstitution;
	
	public static void main(String[] args) throws Exception{
		String smiles = "C(C1C(C(C(C(O1)O)O)O)O)O";		
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer molecule = sp.parseSmiles(smiles);
		HumanMetabolismHelpingFunctions hmhf = new HumanMetabolismHelpingFunctions(0.0, true, 4, "hmdb",true, 30, true);
		ArrayList<Biotransformation> resultArray = hmhf.runCompleteIterations(molecule, 3);
		String outputFileName = "E:/BiotransformerDependencies/Biotransformer3.0/temp_moloutput/testPathBankConnection_result.sdf";
		hmhf.ecb.saveBioTransformationProductsToSdf(resultArray, outputFileName, hmhf.ecb.reactionsHash, false);
		//String smiles_ccooh = "CC(=O)O";
		
		//System.out.println(hmhf.rfpathbank.fcc.retriveFromPathbank(smiles));	
	}
	
	public HumanMetabolismHelpingFunctions(Double scoreThreshold, boolean useDB, int cyp450Mode, String source, boolean checkOverflow, int overflow_threshold, boolean useSubstitution) throws Exception{
		this.useDB = useDB;
		this.useSubstitution = useSubstitution;
		this.cyb = new Cyp450BTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
		this.p2b = new Phase2BTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
		this.ecb	= new ECBasedBTransformer(BioSystemName.HUMAN, useDB, useSubstitution);
		this.hgb = new HGutBTransformer(useDB, useSubstitution);
		//this.mnd_Liver =  new MultivariateDistribution("Liver");
		//this.mnd_ECBased =  new MultivariateDistribution("ECBased");
		//this.mnd_GutSubstrate =  new MultivariateDistribution("GutSubstrate");
		//this.mnd_GutMetabolite =  new MultivariateDistribution("GutMetabolite");
		this.scoreThreshold = scoreThreshold;
		this.cyp450Mode = cyp450Mode;
		if(useDB) {
			this.rfdb = new RetriveFromDB(source, "all", checkOverflow, overflow_threshold, this.useSubstitution);
		}
	}
	
	/**
	 * The user can specify the numIteration he wants to perform
	 * If numIteration = -1, then it means the program will keep running till no more new metabolites will be generated
	 * @param startingCompound
	 * @param numIteration
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> runCompleteIterations(IAtomContainer startingCompound,  int numIteration) throws Exception{	
		System.out.println(this.sg.create(startingCompound));		
		int num_iteration_record = numIteration;
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		Biotransformer bf = new Biotransformer(BioSystemName.HUMAN, this.useDB, this.useSubstitution);	
		//Initialize the substratePool and metabolitePool
		//IAtomContainerSet metabolitesPool = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);		
		IAtomContainerSet substratePool = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		if(this.hp.isPolymer(startingCompound)){
			//substratePool.add(this.hp.convertPolymerToMonomer(startingCompound));
			ArrayList<Biotransformation> monomers_biotransformation = this.hp.converPolymerToMonomer_biotransformation(startingCompound);
			biotransformations.addAll(monomers_biotransformation);
			numIteration = numIteration - 1;
			substratePool.add(getValidStructures(bf.extractProductsFromBiotransformations(biotransformations)));
		}
		else{
			substratePool.addAtomContainer(startingCompound);
		}
		//Create and initialize the substrateProcessed variable to keep track of the processed substrates 
		//If a molecule has been used as substrate and it occurs as metablite in another biotransformation, we should not use it as a substrate in the next iteration
		ArrayList<String> substrateProcessed = new ArrayList<>();	
		IAtomContainerSet novelSubstrateSet;		
		novelSubstrateSet = getNovelSubstrates(substratePool, substrateProcessed);
	
		while(numIteration > 0){
			//Handle substrates that should be retrieved from database
			ArrayList<Biotransformation> currentBiotransformations = new ArrayList<Biotransformation>(); 
			for(int i = 0; i < novelSubstrateSet.getAtomContainerCount(); i++){
				IAtomContainer tempNovelSubstrate = novelSubstrateSet.getAtomContainer(i);
				if(this.useDB && this.rfdb.fcc.retriveFromDB(tempNovelSubstrate)){
					ArrayList<Biotransformation> db_results = this.rfdb.getBiotransformationRetrievedFromDB_setAndstep(tempNovelSubstrate, true);
					System.out.println(this.sg.create(tempNovelSubstrate));
					if(db_results!=null && !db_results.isEmpty()){		
						currentBiotransformations.addAll(db_results);
						novelSubstrateSet.removeAtomContainer(i);
						i = i -1;
					}
				}
			}
			ArrayList<Biotransformation> currentBiotransformations_fromBT = oneCompleteIteration(novelSubstrateSet);			
			currentBiotransformations.addAll(currentBiotransformations_fromBT);
			numIteration--;
			if(!currentBiotransformations.isEmpty()){
				novelSubstrateSet.removeAllAtomContainers();		
				novelSubstrateSet = getValidStructures(bf.extractProductsFromBiotransformations(currentBiotransformations));
				novelSubstrateSet = getNovelSubstrates(novelSubstrateSet, substrateProcessed);		
				substrateProcessed = Utilities.updateProcessedSubstratePool(substrateProcessed, novelSubstrateSet);
				biotransformations.addAll(currentBiotransformations);
			}
			else{
				break;
			}
			
			//Update the metabolitesPool
			//metabolitesPool.add(metabolites_oneIteration);
			//The novel substrates for next iteration should be a subset of the metabolites of the current iteration				
			if(novelSubstrateSet == null || novelSubstrateSet.isEmpty()){
				System.out.println("No more novel substrates at iteration: " + (num_iteration_record - numIteration + 1));
				break;
			}
			System.out.println("num of iterations done: " + (num_iteration_record - numIteration));
			System.out.println("Substrates for next iteration: " + novelSubstrateSet.getAtomContainerCount());
		}
		//Add all biotransformations retrieved from pathbank to the result biotransformation array
		//biotransformations = Utilities.selectByUniqueMetabolites(biotransformations);
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	/**
	 * This function will run EC Based, CYP450, PhaseII and GutMicrobial metabolism parallelly. -- This is called one iteration
	 * The metabolites, except PhaseII metabolites, will then be used as substrates for the next iteration.
	 * @param startingCompounds
	 * @return
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> oneCompleteIteration(IAtomContainerSet startingCompounds) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet currentMetabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<Biotransformation> currentBiots =  new ArrayList<Biotransformation>();
		/**
		 * EC based 
		 */
		IAtomContainerSet ecSubstrates = (IAtomContainerSet) startingCompounds.clone();
		currentBiots = this.ecb.simulateECBasedMetabolismChain(ecSubstrates, true, true, 1, scoreThreshold);
		//currentBiots = this.mnd_ECBased.extractValidBioTransformation_ECBased(currentBiots);
		currentMetabolites.add(this.ecb.extractProductsFromBiotransformations(currentBiots));
		biotransformations.addAll(currentBiots);
		/**
		 * CYP450
		 */	
		IAtomContainerSet cyp450Substrates = (IAtomContainerSet) startingCompounds.clone();		
		currentBiots = this.cyb.predictCyp450BiotransformationsByMode(cyp450Substrates,this.cyp450Mode, true, true, this.scoreThreshold);
		//currentBiots = this.mnd_Liver.extractValidBioTransformation_CYP450(currentBiots);			
		currentMetabolites.add(this.cyb.extractProductsFromBiotransformations(currentBiots));
		biotransformations.addAll(currentBiots);
		/**
		 * Phase II
		 */
		IAtomContainerSet phaseIISubstrates = (IAtomContainerSet) startingCompounds.clone();
		currentBiots = this.p2b.applyPhase2TransformationsChainAndReturnBiotransformations(phaseIISubstrates, true, true, true, 1, this.scoreThreshold);
		IAtomContainerSet p2_metabolites = this.p2b.extractProductsFromBiotransformations(currentBiots);
		//Here we assume that he cunjugated metabolite produced in phase II metabolism is the end product and will not be futher metabolized
		for(int i = 0; i < p2_metabolites.getAtomContainerCount(); i++){
			IAtomContainer one_p2_metabolite = p2_metabolites.getAtomContainer(i);
			one_p2_metabolite.setProperty("isEndProduct", true);
		}
		currentMetabolites.add(p2_metabolites);
		
		biotransformations.addAll(currentBiots);
		/**
		 * Gut Microbial
		 */
		
		IAtomContainerSet gutSubstrates = (IAtomContainerSet) startingCompounds.clone();
		//gutSubstrates = this.mnd_GutSubstrate.isValidGutSubstrates(gutSubstrates);
		currentBiots = this.hgb.applyGutMicrobialMetabolismHydrolysisAndRedoxChain(gutSubstrates, true, true, 1, this.scoreThreshold);
		//currentBiots = this.mnd_GutMetabolite.extractValidBioTransformation_GutMetabolites(currentBiots);
		currentMetabolites.add(this.hgb.extractProductsFromBiotransformations(currentBiots));
		biotransformations.addAll(currentBiots);		
		
		//biotransformations = Utilities.selectByUniqueMetabolites(biotransformations);
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
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
	
	/**
	 * This function will extract unique metabolites from all molecules in the query molecules stored in the IAtomContainerSet
	 * @param allMolecules
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet getUniqueMetabolites(IAtomContainerSet allMolecules, ArrayList<String> processedList) throws Exception{
		IAtomContainerSet cleanedSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		for(int i = 0; i < allMolecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = allMolecules.getAtomContainer(i);
			if(!cyb.isValidMetabolte(oneMole)) continue;
			InChIGenerator inchiGen = inchiFactory.getInChIGenerator(oneMole);
			String inChiKey = inchiGen.getInchiKey();
			if(!processedList.contains(inChiKey)){
				processedList.add(inChiKey);
				cleanedSet.addAtomContainer(oneMole);
			}
		}
		return cleanedSet;
	}
	

	
	/**
	 * This function will explore the substratePool and return the substrates that have not been input into the metabolism iteration
	 * @param substrateSet
	 * @param processedInChiKeyList
	 * @return
	 * @throws Exception
	 */
	public IAtomContainerSet getNovelSubstrates(IAtomContainerSet substrateSet, ArrayList<String> processedInChiKeyList) throws Exception{
		IAtomContainerSet cleanedSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		for(int i = 0; i < substrateSet.getAtomContainerCount(); i++){
			IAtomContainer oneMole = substrateSet.getAtomContainer(i);
			//if(this.isEndProductSMARTS.containEndSubstructure(oneMole)) continue;
			if(this.isEndProductSMARTS.isEndProduct(oneMole)) continue;			
			if(!cyb.isValidMetabolte(oneMole)) continue;
			InChIGenerator inchiGen = inchiFactory.getInChIGenerator(oneMole);
			String inChiKey = inchiGen.getInchiKey();
			if(!processedInChiKeyList.contains(inChiKey)){
				processedInChiKeyList.add(inChiKey);
				cleanedSet.addAtomContainer(oneMole);
			}
		}
		return cleanedSet;
	}
	/**
	 * This function will take the molecules that have at least 5 non-Hydrogen atoms in them
	 * @param allMolecules
	 * @return
	 */
	public IAtomContainerSet getValidStructures(IAtomContainerSet allMolecules) throws Exception{
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);	
		for(int i = 0; i < allMolecules.getAtomContainerCount(); i++){
			if(cyb.isValidMetabolte(allMolecules.getAtomContainer(i))) resultSet.addAtomContainer(allMolecules.getAtomContainer(i));
		}
		return resultSet;
	}
	/**
	 * This function will check if the molecule is a lipid and should gothrough gut microbial metabolism first
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public boolean throughGut(IAtomContainer molecule) throws Exception{
		 ArrayList<ChemicalClassName> chemClasses = ChemicalClassFinder.AssignChemicalClasses(molecule);		 
		 if((chemClasses.contains(ChemicalClassName.ETHER_LIPID) ||
				 chemClasses.contains(ChemicalClassName.GLYCEROLIPID) ||
				 chemClasses.contains(ChemicalClassName.GLYCEROPHOSPHOLIPID) ||
				 chemClasses.contains(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL) ||
				 chemClasses.contains(ChemicalClassName.SPHINGOLIPID) )){
			 return true;
		 }
		 else return false;
	}
	
}
