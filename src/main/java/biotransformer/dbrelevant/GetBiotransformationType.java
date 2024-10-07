package biotransformer.dbrelevant;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;

import biotransformer.transformation.Biotransformation;

public class GetBiotransformationType {
	private final String DB_URL = "jdbc:mysql://162.243.63.130:3306/hmdb";
	private final String USER = "readonly";
	private final String PASS = "tardis";

	public static void main(String[] args) throws Exception{
		GetBiotransformationType gbt = new GetBiotransformationType();
		try(
				Connection conn = DriverManager.getConnection(gbt.DB_URL, gbt.USER, gbt.PASS);
				Statement stmt = conn.createStatement();
			){		      
				System.out.println(gbt.isCYP450Biotransformation(stmt, 1439));
				System.out.println(gbt.isPhaseIIBiotransformation(stmt, 2348));
				System.out.println(gbt.isMicrobialBiotransformation(stmt, 2348));
			 } catch (Exception e) {
					e.printStackTrace();
				}
	}
	public boolean isQueriedBiotransformation(Statement stmt, Integer reactionID, String query_biotransformation_type) throws Exception{
		if(query_biotransformation_type.equals("all")) {
			return true;
		}
		else if(query_biotransformation_type.equals("cyp450")) {
			return isCYP450Biotransformation(stmt, reactionID);
		}
		else if(query_biotransformation_type.equals("phaseII")) {
			return isPhaseIIBiotransformation(stmt, reactionID);
		}
		else if(query_biotransformation_type.equals("microbial")) {
			return isMicrobialBiotransformation(stmt, reactionID);
		}
		//Dr. Wishart suggested that the microbial metabolism should cover almost all the human metabolism and gut microbial metabolism reactions
		else if(query_biotransformation_type.equals("envmicrobial")) {
			return true;
			
		}
		else throw new Exception("Please specify the biotransformation type you want to query when calling retrieving from db feature. The types are all, cyp450, phaseII and microbial");
	}
	/**
	 * This function checks if the retrieved biotransformation is a CYP450 biotransformation or not by:
	 * Get the enzyme A that catalyzes the biotransformation, then get the gene_name of enzyme A.
	 * If the gene_name contains keyword "CYP", then it's a CYP450 enzyme and thus the biotransformation is catalyzed by CYP450 enzymes.
	 * @param stmt
	 * @param reactionID
	 * @return
	 * @throws Exception
	 */
	public boolean isCYP450Biotransformation(Statement stmt, Integer reactionID) throws Exception{
		String query = "select gene_name, protein_id from reaction_enzymes, proteins where reaction_id = '%d' and proteins.id = protein_id";
		query = String.format(query, reactionID);
		ResultSet rs = stmt.executeQuery(query);
		while(rs.next()) {
			//Integer protein_id = rs.getInt("protein_id");
			String gene_name = rs.getString("gene_name");
			if(gene_name == null) return false;
			if(gene_name.contains("CYP") || gene_name.contains("cyp") || gene_name.contains("Cyp")) {
				rs.close();
				return true;
			}
		}
		rs.close();
		return false;
	}

	/**
	 * This function checks if the retrieved biotransformation is a PhaseII biotransformatin or not by:
	 * Get the name of the enzyme A that catalyze the biotransformation specified by the input reactionID.
	 * If the name contains any of the following keywords:  sulfotransferase, methyltransferase, N-acetyltransferase, UDP-glucuronosyltransferase, 	glutathione S-transferase, Glycine N-acyltransferase
	 * then it is a biotransformation catalyzed by PhaseII enzyme
	 * @param stmt
	 * @param reactionID
	 * @return
	 * @throws Exception
	 */
	public boolean isPhaseIIBiotransformation(Statement stmt, Integer reactionID) throws Exception{
		String query = "select name, protein_id from reaction_enzymes, proteins where reaction_id = '%d' and proteins.id = protein_id";
		query = String.format(query, reactionID);
		ResultSet rs = stmt.executeQuery(query);
		while(rs.next()) {
			Integer protein_id = rs.getInt("protein_id");
			String enzyme_name = rs.getString("name");
			if(isValidPhaseIIEnzyme(enzyme_name)) {
				rs.close();
				return true;
			}
		}		
		rs.close();
		return false;
	}
	
	/**
	 * sulfotransferase, methyltransferase, N-acetyltransferase, UDP-glucuronosyltransferase, 	glutathione S-transferase, Glycine N-acyltransferase
	 * @param enzyme
	 * @return
	 */
	public boolean isValidPhaseIIEnzyme(String enzyme) {
		if(enzyme.contains("sulfotransferase") || enzyme.contains("methyltransferase") || enzyme.contains("N-acetyltransferase") ||
		   enzyme.contains("UDP-glucuronosyltransferase") || enzyme.contains("glutathione S-transferase") || enzyme.contains("Glycine N-acyltransferase")) {
			return true;
		}
			
		return false;
	}
	
	public boolean isMicrobialBiotransformation(Statement stmt, Integer reactionID) throws Exception{
		String query = "select name, organism,  protein_id from reaction_enzymes, proteins where reaction_id = '%d' and proteins.id = protein_id";
		query = String.format(query, reactionID);
		ResultSet rs = stmt.executeQuery(query);
		while(rs.next()) {
			String organism = rs.getString("organism");
			String enzyme_name = rs.getString("name");
			//Treat this case as true for now. This should be rare.
			if(organism == null) {
				//System.out.println(reactionID + ", " + organism + ", " + enzyme_name);
				return true;
			}
			if(organism.contains("Homo") || organism.contains("homo") || organism.contains("HOMO")) {
				rs.close();
				return false;
			}
		}		
		rs.close();
		return true;
	}
	public boolean isMicrobialEnzyme(String enzyme) {
		return false;
	}
	
}
