import java.util.HashMap;

public class PloidyManager {
	public HashMap <String,Integer> ploidies;
	
	public PloidyManager(){
		ploidies = new HashMap<String,Integer>();
	}
	
	public Integer put( String contig,int ploidy){

		return ploidies.put(contig, ploidy);
	}
	
	public Integer get(String contig){
		return ploidies.get(contig);
	}
	
}
