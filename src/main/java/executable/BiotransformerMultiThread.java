package executable;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

//import com.opencsv.CSVReader;
//import com.opencsv.CSVWriter;
import java.io.Reader;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.csv.CSVPrinter;

import biotransformer.btransformers.SimulateHumanMetabolism;
import biotransformer.utils.HumanSuperBioTransformer;

public class BiotransformerMultiThread {
    /**
     * Main class to run multi-thread
     * numCores, prediction_mode (useless now), inputfile, time_max, cyp450Mode, iteration, useDB
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception{


        IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
        SmilesParser	smiParser		= new SmilesParser(builder);


        // List<Future<Boolean>> list = new ArrayList<Future<Boolean>>();
        HashMap<String,Future<Boolean>> maps = new HashMap<String,Future<Boolean>>();

        int numberCore = Integer.valueOf(args[0]);
        String mode = args[1];
        if(mode == null) {
            System.out.println("Set the mode (superbio;allHuman)");
            System.exit(0);
        }

        String inputfileName = args[2];
        Integer exception_time = Integer.valueOf(args[3]);
        System.out.println(String.format("Current exception time is %d", exception_time));

        String current_dir = System.getProperty("user.dir");

        System.out.println("number of core: " + Integer.toString(numberCore));
        ExecutorService executor = Executors.newFixedThreadPool(numberCore);
        int current_file_index = Integer.valueOf(inputfileName.split("_")[1].replace(".csv",""));

        String inputFile = String.format("%s/%s/%s", current_dir,"inputFolder",inputfileName);
        String outputFile = String.format("%s/%s_%d.csv", current_dir,"inputFolder/BioTransformerExeption",current_file_index);
        Reader csvReader = new FileReader(inputFile);
        Iterable<CSVRecord> records = CSVFormat.EXCEL.parse(csvReader);

        CSVPrinter csvFilePrinter = null;
        CSVFormat csvFileFormat = CSVFormat.EXCEL.withHeader();
        FileWriter fileWriter = new FileWriter(outputFile);
        CSVPrinter StatusWriter = new CSVPrinter(fileWriter, csvFileFormat);

        //CSVReader csvReader = new CSVReader(new FileReader(String.format("%s/%s/%s", current_dir,"inputFolder",inputfileName)));
        //CSVWriter StatusWriter = new CSVWriter(new FileWriter(String.format("%s/%s_%d.csv", current_dir,"inputFolder/BioTransformerExeption",current_file_index)));
        String[] nextLine;
        // csvReader.readNext(); // skip header
        int cyp450mode = 3;//We chose mode 3: both cyProduct and original biotransformer as default value
        //while ((nextLine = csvReader.readNext()) != null) {
        for (CSVRecord record : records) {

            if (record != null) {
                //String molID = nextLine[0];
                String molID = record.get(0);
                //String structure = nextLine[2];
                String structure = record.get(2);
                String smiles = structure;

                String fileName = String.format("%s/temp_moloutput/molouput_%s.sdf", current_dir,molID);

                IAtomContainer singleInput = smiParser.parseSmiles(smiles);
                InChIGeneratorFactory inchiFa = InChIGeneratorFactory.getInstance();
                System.out.println(smiles);
                System.out.println(inchiFa.getInChIGenerator(singleInput).getInchiKey());
                InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
                String inChikey = inchiFactory.getInChIGenerator(singleInput).getInchiKey();
                System.out.println(inChikey);

                Future<Boolean> future = executor.submit(new CallableSimulateHuman(singleInput,fileName, cyp450mode, 1, false));
                maps.put(molID, future);
            }

        }


        for(String key : maps.keySet()) {
            try {
                // https://docs.oracle.com/javase/6/docs/api/java/util/concurrent/Future.html#get%28long,%20java.util.concurrent.TimeUnit%29
                //Future<Boolean> f = maps.get(key);
                //if(f.isDone()){
                boolean result = maps.get(key).get(exception_time, TimeUnit.SECONDS);
                if(!result) {
                    StatusWriter.printRecord(new String[] {key,"ProgramException"});
                }
                //}

            }catch (TimeoutException e) {
                maps.get(key).cancel(true);
                StatusWriter.printRecord(new String[] {key,"TimeoutException"});

            }catch (InterruptedException e) {
                StatusWriter.printRecord(new String[] {key,"InterruptedException"});
            } catch (ExecutionException e) {
                e.printStackTrace();
                StatusWriter.printRecord(new String[] {key,"ExecutionException"});
            }catch (Exception e) {
                e.printStackTrace();
                StatusWriter.printRecord(new String[] {key,"Exception"});
            }
        }


        executor.shutdownNow();
        csvReader.close();
        StatusWriter.close();
    }
    public static void runMultiThread(String[] args) throws Exception{
    	IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
        SmilesParser	smiParser		= new SmilesParser(builder);
        // List<Future<Boolean>> list = new ArrayList<Future<Boolean>>();
        HashMap<String,Future<Boolean>> maps = new HashMap<String,Future<Boolean>>();

        int numberCore = Integer.valueOf(args[0]);

        String inputfileName = args[1];
        Integer exception_time = Integer.valueOf(args[2]);
        System.out.println(String.format("Current exception time is %d", exception_time));

        String current_dir = System.getProperty("user.dir");

        System.out.println("number of core: " + Integer.toString(numberCore));
        ExecutorService executor = Executors.newFixedThreadPool(numberCore);
        int current_file_index = Integer.valueOf(inputfileName.split("_")[1].replace(".csv",""));

        String inputFile = String.format("%s/%s/%s", current_dir,"inputFolder",inputfileName);
        String outputFile = String.format("%s/%s_%d.csv", current_dir,"inputFolder/BioTransformerExeption",current_file_index);
        Reader csvReader = new FileReader(inputFile);
        Iterable<CSVRecord> records = CSVFormat.EXCEL.parse(csvReader);

        CSVPrinter csvFilePrinter = null;
        CSVFormat csvFileFormat = CSVFormat.EXCEL.withHeader();
        FileWriter fileWriter = new FileWriter(outputFile);
        CSVPrinter StatusWriter = new CSVPrinter(fileWriter, csvFileFormat);

        //CSVReader csvReader = new CSVReader(new FileReader(String.format("%s/%s/%s", current_dir,"inputFolder",inputfileName)));
        //CSVWriter StatusWriter = new CSVWriter(new FileWriter(String.format("%s/%s_%d.csv", current_dir,"inputFolder/BioTransformerExeption",current_file_index)));
        String[] nextLine;
        // csvReader.readNext(); // skip header
        int cyp450mode = Integer.parseInt(args[3]);//We chose mode 3: both cyProduct and original biotransformer as default value
        int iteration = Integer.parseInt(args[4]);
        boolean useDB = Boolean.parseBoolean(args[5]);
        //while ((nextLine = csvReader.readNext()) != null) {
        for (CSVRecord record : records) {
            if (record != null) {
                //String molID = nextLine[0];
                //String structure = nextLine[2];
                //String smiles = structure;
                String molID = record.get(0);
                String structure = record.get(2);
                String smiles = structure;
            	try{	            
	                String folderName = String.format("%s/temp_moloutput/%s/", current_dir,current_file_index);
	                File folder =  new  File(folderName);
	                if(!folder.exists() || !folder.isDirectory()){
	                	folder.mkdir();
	                }
	                String fileName = String.format("%s/molouput_%s.sdf", folderName,molID);
	
	                IAtomContainer singleInput = smiParser.parseSmiles(smiles);
	                InChIGeneratorFactory inchiFa = InChIGeneratorFactory.getInstance();
	                //System.out.println(smiles);
	                //System.out.println(inchiFa.getInChIGenerator(singleInput).getInchiKey());
	                InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
	            	String inChikey = inchiFactory.getInChIGenerator(singleInput).getInchiKey();
	            	//System.out.println(inChikey);
	                
	                Future<Boolean> future = executor.submit(new CallableSimulateHuman(singleInput,fileName, cyp450mode, iteration, useDB));
	                maps.put(molID, future);
            	}catch (Exception e){
            		System.out.println("Failed processing molecule: " + smiles);
            		continue;
            	}
            }

        }


        for(String key : maps.keySet()) {
            try {
                // https://docs.oracle.com/javase/6/docs/api/java/util/concurrent/Future.html#get%28long,%20java.util.concurrent.TimeUnit%29
            	//Future<Boolean> f = maps.get(key);
            	//if(f.isDone()){
            		boolean result = maps.get(key).get(exception_time, TimeUnit.SECONDS);
            		if(!result) { 
            			StatusWriter.printRecord(new String[] {key,"ProgramException"});
            			}
            	//}

            }catch (TimeoutException e) {
                maps.get(key).cancel(true);
                StatusWriter.printRecord(new String[] {key,"TimeoutException"});

            }catch (InterruptedException e) {
                StatusWriter.printRecord(new String[] {key,"InterruptedException"});
            } catch (ExecutionException e) {
                e.printStackTrace();
                StatusWriter.printRecord(new String[] {key,"ExecutionException"});
            }catch (Exception e) {
                e.printStackTrace();
                StatusWriter.printRecord(new String[] {key,"Exception"});
            }
        }


        executor.shutdownNow();
        csvReader.close();
        StatusWriter.close();
    }
    /**
     * Callable for AllHuman mode
     */
    private static class CallableSimulateHuman implements Callable<Boolean> {
        private IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
        private IAtomContainer mole = builder.newAtomContainer();
        private String fileName = new String();
        private Integer modeNumber = new Integer(0);
        private int iteration = 0;
        private boolean useDB = false;
        public CallableSimulateHuman(IAtomContainer mole, String fileName, Integer modeNumber, int iteration, boolean useDB)
        {
            this.mole = mole;
            this.fileName = fileName;
            this.modeNumber = modeNumber;
            this.iteration = iteration;
            this.useDB = false;
        }

        public Boolean call() throws Exception{
            try {
            	String thread = Thread.currentThread().getName();  
            	System.out.println("Thread: " + thread + " is running");
            	//HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();
                //hsbt.simulateAllHumanMetabolismAndSavetoSDF(this.mole, this.fileName, 0.5, false, this.modeNumber);
                SimulateHumanMetabolism shm = new SimulateHumanMetabolism(this.modeNumber, this.useDB, "hmdb", false, 30, false);
                shm.simulateHumanMetabolismAndSaveToSDF(this.mole, this.fileName, iteration, false);
                System.out.println("Thread: " + thread + " finished");
            	return true;
            }
            catch(Exception e) {
            	e.printStackTrace();
                return false;
            }

        }
    }
    
    
    /**
     * Callable for AllHuman mode
     */
    private static class CallableAllHuman implements Callable<Boolean> {
        private IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
        private IAtomContainer mole = builder.newAtomContainer();
        private String fileName = new String();
        private Integer modeNumber = new Integer(0);
        private boolean useDB = false;
        private boolean useSubstitution = false;
        public CallableAllHuman(IAtomContainer mole, String fileName, Integer modeNumber)
        {
            this.mole = mole;
            this.fileName = fileName;
            this.modeNumber = modeNumber;

        }

        public Boolean call() throws Exception{
            try {
                HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer(useDB, this.useSubstitution);
                hsbt.simulateAllHumanMetabolismAndSavetoSDF(this.mole, this.fileName, 0.5, false, this.modeNumber);
                return true;
            }
            catch(Exception e) {
                return false;
            }

        }
    }


    /**
     * Callable for Superbio mode
     */
    private static class CallableSuperbio implements Callable<Boolean> {
        private IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
        private IAtomContainer mole = builder.newAtomContainer();
        private String fileName = new String();
        private boolean useDB = false;
        private boolean useSubstitution = false;
        public CallableSuperbio(IAtomContainer mole, String fileName)
        {
            this.mole = mole;
            this.fileName = fileName;

        }

        public Boolean call() throws Exception{
            try {
                HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer(useDB, useSubstitution);
                hsbt.simulateHumanSuperbioMetabolismAndSaveToSDF(this.mole, this.fileName, false);
                return true;
            }
            catch(Exception e) {
                return false;
            }



        }
    }
}
