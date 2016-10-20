import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.stream.Stream;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class PloidyPhaser {

	static String ploidyFile = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//PhasingAndVCFs//PloidyEstimationTest.txt";
	static String vcfFolder = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//Simulated//PhasingAndVCFs//MOVECHROM1OUT";
	static List<File> bamPaths = new ArrayList<File>();
	static HashMap <String,Integer> readIndexes = new HashMap<String,Integer>(16336669);//index of reads
	static HashMap ploidies = new HashMap();
	static List<Path> vcfPaths = new ArrayList<Path>();
	
	public static void main(String[] args) {

		getFiles();
		for (int i = 0; i < vcfPaths.size(); i++) {//for each contig, a vcf file is associated
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
			bamPaths.add(new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//ConPADE-master//Simulated_input_bams//sorted_200bp.bam"));
			bamPaths.add(new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//ConPADE-master//Simulated_input_bams//sorted_500bp.bam"));
			//getReadIndexes();
			getChr1ReadIndexes();
		} catch (Exception e) {
			System.out.println("Error reading ploidy estimation file");
			e.printStackTrace();

		}

	}

	
	private static void getReadIndexes() throws FileNotFoundException, UnsupportedEncodingException {

		int readInd=1;//read index
		SAMFileReader inputSam;
		String libraryName;
		int ct=0;
		for (int bi=0;bi<PloidyPhaser.bamPaths.size();bi++){
			
			inputSam = new SAMFileReader(PloidyPhaser.bamPaths.get(bi));
			
			//get fileNameHeader
			libraryName="."+bamPaths.get(bi).getName().substring(0,bamPaths.get(bi).getName().lastIndexOf('.'));
			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter = inputSam.iterator();
			String chromName = "";

			while (iter.hasNext()  ) {//iterates the sam file
				SAMRecord rec = iter.next();
				chromName=rec.getReadName()+libraryName;//add library to chromosome name (makes it unique over all libraries)
				//for the moment keep only reads from chromosome 1
				if(rec.getReadName().substring(20, 23).equals("_1_") && !readIndexes.containsKey(chromName) ){
					readIndexes.put(chromName, readInd++);	
				}//else if (!rec.getReadName().substring(20, 23).equals("_1_"))	break;	//ct++;		
			}
			iter.close();
			inputSam.close();
			System.out.println("ct "+(ct)+" libraryName:"+libraryName+" nbReads:"+readIndexes.size());
			
		}
		
		
		/*
		PrintWriter writer = new PrintWriter("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//ConPADE-master//Simulated_input_bams//"+"chrom1ReadsIndexes.txt", "UTF-8");
		
		for (HashMap.Entry<String, Integer> entry : readIndexes.entrySet()) {
			writer.println(entry.getKey() + "\t" + entry.getValue());
		}
		writer.close();
		*/
	}
	
	//method to spped up coding- REMOVE WHEN FINISHED and replace with CORRECTED getReadIndexes() (corrected means remove all references to chromosome 1)
	private static void getChr1ReadIndexes()  {

		try {
			Scanner sc = new Scanner(new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//ConPADE-master//Simulated_input_bams//chrom1ReadsIndexes.txt"));
			
			while (sc.hasNextLine()){
				readIndexes.put(sc.next(), sc.nextInt());
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	}
	
	
	
}
