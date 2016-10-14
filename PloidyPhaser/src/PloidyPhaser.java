import java.io.File;
import java.io.FileNotFoundException;
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

	static String ploidyFile = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//PhasingAndVCFs//PloidyEstimationTest.txt";
	static String vcfFolder = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//PhasingAndVCFs//MOVECHROM1OUT";
	static List<File> bamPaths = new ArrayList<File>();

	static HashMap ploidies = new HashMap();
	static List<Path> vcfPaths = new ArrayList<Path>();

	public static void main(String[] args) {

		getFiles();
		for (int i = 0; i < vcfPaths.size(); i++) {
			try {
				VCFManager vcfMan = new VCFManager(vcfPaths.get(i).toFile());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

	}

	private static void getFiles() {
		String line = "";
		// System.out.println("Paths.get(vcfFolder).toString():"+Paths.get(vcfFolder).toString()
		// );

		try {
			// get ploidies
			Scanner sc = new Scanner(new File(ploidyFile));
			while (sc.hasNextLine()) {
				// n++;
				line = sc.nextLine();
				String[] parts = line.split("\\s+");
				ploidies.put(parts[0], Integer.parseInt(parts[1]));// key:chromosome
																	// name/value:ploidy
			}

			// get vcf Files
			Stream<Path> paths = Files.walk(Paths.get(vcfFolder));
			paths.forEach(filePath -> {
				// consider only file in input folder (ignore subfolders)
				if (Files.isRegularFile(filePath)
						&& filePath.toString().substring(filePath.toString().lastIndexOf(".") + 1).equals("vcf")
						&& filePath.getParent().toString().equals(Paths.get(vcfFolder).toString())) {
					vcfPaths.add(filePath);
					// System.out.println("filepath:"+filePath );
				}
			});
			//get bams
			bamPaths.add(new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//BAMs//200bp.sam"));
		} catch (Exception e) {
			System.out.println("Error reading ploidy estimation file");

		}

	}

}
