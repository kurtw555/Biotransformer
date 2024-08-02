/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.biosystems;

import java.io.FileNotFoundException;
import java.io.IOException;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import exception.BioTransformerException;

public class Human extends BioSystem {

	public Human(ObjectMapper mapper) throws JsonParseException, JsonMappingException, 
	FileNotFoundException, IOException, BioTransformerException {
		super(BioSystemName.HUMAN, mapper);
	}

}
