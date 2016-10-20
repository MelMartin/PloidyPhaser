
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;


public class VCFManager {
	
	String header="";//vcf Header
	List<VariationData> cleanLines = new ArrayList<VariationData>();
	String outputCleanFile="" ;
	File vcfFile;
	VariationsManager varMan=new VariationsManager(this);
	
	
	
	
	//----------------------- AM I GOING TO NEED THESE????--------------------------
	
	List<Integer> refInsertsPos=new ArrayList<Integer>();//keeps track of inserts pos into the reference allele (each time 'alt' is longer than 'ref')
	List<Integer> nbInsertsPos=new ArrayList<Integer>();//keeps track of nb of inserts into the reference allele (each time 'alt' is longer than 'ref')
	
	//-----------------------ERASE IF NOT USED -------------------------------------


	public VCFManager(File vf) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		
		vcfFile=vf;
		
		vcfVariantShortExtractor(vcfFile.getAbsolutePath());//extracts only the repeated variants in a short format: chromName, pos, ref, alt
		varMan.constructVariantMatrix();
		//varMan.fillVariantMatrix();

	}

	


	public void vcfVariantShortExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

		String line="";
		
		Scanner sc = new Scanner(new File(inputFile));

		int ct = 0;
		String chrom="";
		int pos = 0;
		int svLength;// length of the sv
		int end;
		Boolean isPloidySolved=false;
		Boolean isVariation=false;
		int ploidy=0;
		String ref = "";// reference allele
		String alt=".";// alternative allele
		List<String> infoFields=new ArrayList();

		String format="GT";
		String sample="";
		VariationData currentVariationData=new VariationData() ;

		//String next="";
		try {
			// get the header		
			line=sc.nextLine();
			while (line.substring(0, 1).equals("#")){	
				if(!line.substring(0,6).equals("#CHROM")){
					header+=line ;	
				}else header+="#CHROM					POS	REF	ALT	SAMPLE";//but adapt the fields to the new short format
				line=sc.nextLine();				
				if(line.substring(0, 1).equals("#"))header+="\n";
				ct++;
			}

			//reset scaner after the header;
			sc = new Scanner(new File(inputFile));
			for(int l=0;l<ct;l++){
				sc.nextLine();//skip headers line
			}
			ct=0;
			

			//GET THE VCF VALUES
			chrom=sc.next();
			//solve ploidy (general for all lines)
			if (!isPloidySolved ){
				ploidy=(int) PloidyPhaser.ploidies.get(chrom);
				isPloidySolved=true;
			}
			//solve sample (general for all lines)
			if (ploidy==1)sample="1";
			else sample="0/1";
			if( ploidy>2){					
				for (int p=2;p<ploidy;p++){
					sample+="/0";//just for the right number of copies, it doesn't represent the real phasing positions
				}
			}
			//get each specific line
			while (sc.hasNextLine()  && pos<10010 /* */) {//test with only the first 10.000 positions

				pos = Integer.parseInt(sc.next());	// get pos	
				sc.next();// skip id
				ref = sc.next();// get reference allele
				alt = sc.next();// get alternative allele
				sc.next();//skip qual		
				sc.next();//skip nextFilter
				infoFields = Arrays.asList(sc.next().split(";"));// split all fields in the info line				
				format=sc.next();
				sc.next();//skip sample
	
				//buid the line that is going to be written
				if (!infoFields.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )

					svLength = Integer.parseInt(infoFields.get(1).substring(6, infoFields.get(1).length())) ;
					end=Integer.parseInt(infoFields.get(2).substring(4, infoFields.get(2).length())) ;
					System.out.println("pos:"+pos+" SV infoFields:"+infoFields+ " LENGTH:" + svLength+" FORMAT:"+format);
					currentVariationData = new VariationData(chrom, pos ,ref,alt,sample);
					isVariation=true;
					//svLength = 0;
				} else {// ...regular vcf line
					currentVariationData =new VariationData(chrom, pos ,ref,alt,sample);		
				}
				

				//write out the builded lines
				if((isVariation || !alt.equals("."))){
					cleanLines.add(currentVariationData);
					isVariation=false;
					varMan.registerVariation(currentVariationData);
					if (alt.length()>ref.length()){
						refInsertsPos.add(pos);
						nbInsertsPos.add(alt.length()-ref.length());
					}
				}

				if (sc.hasNextLine() ) { // 'if' to avoid error at end of file
					currentVariationData =  new VariationData();
					chrom=sc.next();// contig name 
				}else {
					cleanLines.add(currentVariationData);//adds last line
				}
				ct++;		
			}

			printCleanLines("Short");

			if (sc != null)
				sc.close();
			
			//CreateVariationVertex();//creates the vertex that are going to be used to define the phasing Matrix 'varMat'

		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentVariationData.outString());
		}

	}



	public void printCleanLines(String apendix) throws FileNotFoundException {

		System.out.println("outputCleanFile:"+vcfFile.getParent()+"\\Clean\\"+vcfFile.getName().substring(0, vcfFile.getName().lastIndexOf('.'))+apendix+".vcf");
		outputCleanFile = vcfFile.getParent()+"\\Clean\\"+vcfFile.getName().substring(0, vcfFile.getName().lastIndexOf('.'));

		outputCleanFile +=apendix+".vcf";

		PrintStream stdout = System.out;
		PrintStream myConsole = null;

		try {
			myConsole = new PrintStream(new File(outputCleanFile));
			System.setOut(myConsole);
			System.out.println(header);
			for (int i = 0; i < (cleanLines.size()-1); i++) {
				System.out.println(cleanLines.get(i).outString());
			}
			System.out.print(cleanLines.get(cleanLines.size()-1).outString());//last line without endline
			myConsole.close();
		} catch (Exception e) {
			System.err.println("Error trying to write outputHumanFile");
			e.printStackTrace();
		}  finally {
			if (myConsole != null) {	        	
				myConsole.close();
				System.setOut(stdout);  
			}
		}
		System.setOut(stdout);  
	}
	
	/*
	private void CreateVariationVertex() {
		VariationData var;
		
		for (HashMap.Entry<VariationData,Integer> entry : varMan.variations.entrySet()) {
			var=entry.getKey();
			if (var.ref.length()>1){
				var.printOut();
				//int dif=var.ref.length()
			}
		}
		//TODO with real size
		//varMat=new int[variationsREAL.size()][variationsREAL.size()];
		
	}
*/
	
}
