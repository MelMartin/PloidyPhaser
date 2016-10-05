import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Stream;


public class PloidyPhaser {

	static String ploidyFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//PhasingAndVCFs//PloidyEstimationTest.txt";
	static String vcfFolder ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//PhasingAndVCFs//MOVECHROM1OUT";

	static HashMap  ploidies= new HashMap();
	static List<Path> vcfPaths=new ArrayList<Path>();
	
	public static void main(String[] args) {
		String line="";
		try{
			Scanner sc = new Scanner(new File(ploidyFile));			
			while (sc.hasNextLine()){
				//n++;
				line = sc.nextLine();
				String[] parts=line.split("\\s+");
				ploidies.put(parts[0],Integer.parseInt(parts[1]));
			}
			
			Stream<Path> paths = Files.walk(Paths.get(vcfFolder));

			    paths.forEach(filePath -> {
			        if (Files.isRegularFile(filePath) && filePath.toString().substring(filePath.toString().lastIndexOf(".") + 1).equals("vcf")) {
			        	vcfPaths.add(filePath);
			        	System.out.println(filePath );
			        }
			    });
			for (int i=0;i<vcfPaths.size();i++){
				VCFManager vcfMan=new VCFManager(vcfPaths.get(i).toFile());
			}
			
		}catch (Exception e) {
			System.out.println("Error reading ploidy estimation file" );

		}

	}

	

}
