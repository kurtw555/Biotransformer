package executable;

import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class TestMultiThread {
	@SuppressWarnings("rawtypes")
	public static void main(String[] args) throws Exception{
		HashMap<String,Future<Boolean>> maps = new HashMap<String,Future<Boolean>>();
		ExecutorService executor = Executors.newFixedThreadPool(3);
		Future<Boolean> future_1 = executor.submit(new CallableSimulateHuman(10));
		Future<Boolean> future_2 = executor.submit(new CallableSimulateHuman(5));
		maps.put("2", future_2);
		maps.put("1", future_1);
		//System.out.println(future_2.get());
		for(String key : maps.keySet()){
			boolean result = maps.get(key).get();
			System.out.println(key + ": " + result);
		}
		executor.shutdown();
	}
	
    /**
     * Callable for AllHuman mode
     */
    private static class CallableSimulateHuman implements Callable<Boolean> {
    	int counter = 10;
        public CallableSimulateHuman(int counter){
           this.counter = counter;
        }

        public Boolean call() throws Exception{
        	String thread = Thread.currentThread().getName();           	
        	for(int i = 0; i < counter; i++){
        		System.out.println(thread + ": " + i);
        		Thread.sleep(1000);
        	}
        	return true;
        }
    }
}
