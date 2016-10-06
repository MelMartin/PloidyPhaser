
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;


public class VCFManager {
	List<String> humanLines = new ArrayList<String>();
	String outputCleanFile="" ;


	
	public VCFManager(File vcfFile) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		//solve the paths
		System.out.println("outputHumanFile:"+vcfFile.getParent()+"\\"+vcfFile.getName().substring(0, vcfFile.getName().lastIndexOf('.'))+"Clean.vcf");
		outputCleanFile = vcfFile.getParent()+"\\"+vcfFile.getName().substring(0, vcfFile.getName().lastIndexOf('.'))+"Clean.vcf";
		
		vcfExtractor(vcfFile.getAbsolutePath());	
	}
	
	public void vcfExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

		
		//outputMatrixFile = outputFileRoot  + "\\"+ endFix + "MatFileOut.txt";
		//System.out.println("outputHumanFile" +outputHumanFile + "   outputMatrixFile" + outputMatrixFile );		
		// We are interested in the 6th and 7th column:
		// 6th="Count of As, Cs, Gs, Ts at locus"
		// 7th="="Percentage of As, Cs, Gs, Ts weighted by Q & MQ at locus"
		String line = "";
		Scanner sc = new Scanner(new File(inputFile));
		
		int ct = 0;
		String chrom="";
		int prevPos=0;
		int pos = 0;
		int depth = 0;// depth
		int bcA = 0;
		int bcC = 0;
		int bcG = 0;
		int bcT = 0;
		int qpA = 0;
		int qpC = 0;
		int qpG = 0;
		int qpT = 0;
		int sv;// structural variance (boolean)
		int svLength;// length of the sv
		int end;
		Boolean isPloidySolved=false;
		int ploidy=0;
		String ref = "";// reference allele
		String alt=".";// alternative allele
		List<String> infoFields;
		List<String> baseCalls;
		List<String> qPercentages;
		String format="GT";
		String sample="";
		String nextFilter;
		String currentHumanLine = "";
		String previousHumanLine = "";

		String next;
		try {
			// skip the header			
			do{
				line = sc.nextLine();
				System.out.println(line);
				next=sc.next();
			}while (next.substring(0, 1).equals("#")); 
			chrom=next;
			
			//get the values
			while (sc.hasNextLine() /* && ct < 34000 */) {
				next=sc.next();				
				pos = Integer.parseInt(next);// get pos
				sc.next();// skip id
				ref = sc.next();// get reference allele
				alt = sc.next();// get alternative allele
				sc.next();// skip qual
				nextFilter = sc.next();// get the filter call. Ambiguous and
				// Deletions are the interesting ones

				infoFields = Arrays.asList(sc.next().split(";"));// split all
				// fields in
				// the info
				// line
				format=sc.next();
				
				if (!isPloidySolved ){
					ploidy=(int) PloidyPhaser.ploidies.get(chrom);
					isPloidySolved=true;
				}
				sample=sc.next();
				
				if( ploidy>2){
					for (int p=2;p<ploidy;p++){
						sample+="/0";
					}
				}
				
				if (!infoFields.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )
					//System.out.println("SV");
					depth = 0;
					bcA = 0;
					bcC = 0;
					bcG = 0;
					bcT = 0;
					qpA = 0;
					qpC = 0;
					qpG = 0;
					qpT = 0;
					sv = 1;
					svLength = Integer.parseInt(infoFields.get(1).substring(6, infoFields.get(1).length())) ;
					end=Integer.parseInt(infoFields.get(2).substring(4, infoFields.get(2).length())) ;
					
					System.out.println("pos:"+pos+" SV infoFields:"+infoFields+ " LENGTH:" + svLength+" FORMAT:"+format);
					currentHumanLine = chrom+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" + nextFilter + "\tSVTYPE="
							+ infoFields.get(0) + ";SVLEN=" + svLength+";END="+end+"\t"+format+"\t"+sample;

				} else {// ...regular vcf line
					depth = Integer.parseInt(infoFields.get(0).substring(3, infoFields.get(0).length()));
					baseCalls = Arrays
							.asList(infoFields.get(5).substring(3, infoFields.get(5).length()).split(","));

					bcA = Integer.parseInt(baseCalls.get(0));
					bcC = Integer.parseInt(baseCalls.get(1));
					bcG = Integer.parseInt(baseCalls.get(2));
					bcT = Integer.parseInt(baseCalls.get(3));

					qPercentages = Arrays
							.asList(infoFields.get(6).substring(3, infoFields.get(6).length()).split(","));

					qpA = Integer.parseInt(qPercentages.get(0));
					qpC = Integer.parseInt(qPercentages.get(1));
					qpG = Integer.parseInt(qPercentages.get(2));
					qpT = Integer.parseInt(qPercentages.get(3));
					sv = 0;
					svLength = 0;
					currentHumanLine = chrom+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" + nextFilter + "\t";
					for (int f=0;f<(infoFields.size()-2);f++){
						currentHumanLine += infoFields.get(f)+";";
					}
					currentHumanLine +=   infoFields.get(infoFields.size()-1)+"\t"+format+ "\t"+sample;
				}
				if(pos!=prevPos && prevPos!=0){
					humanLines.add(previousHumanLine);
				}
				prevPos=pos;
				previousHumanLine = currentHumanLine;

				line = sc.nextLine();
				if (sc.hasNextLine()) { // 'if' to avoid error at end of file
					chrom=sc.next();// contig name 
				}
				ct++;
			}

			printHumanLines();

			if (sc != null)
				sc.close();
		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentHumanLine);
		}

	}

	public int bitSequence(String s) {
		char[] cseq = new char[s.length()];
		int chASCII;
		for (int i = 0; i < s.length(); i++) {
			chASCII = (int) Character.toUpperCase(s.charAt(i));
			// System.out.println("processing:"+chASCII+":"+s.charAt(i)+":");
			switch (chASCII) {
			case 65://'A'
				cseq[i] = '1';
				break;
			case 67://'C'
				cseq[i] = '2';
				break;
			case 71://'G'
				cseq[i] = '3';
				break;
			case 84://'T'
				cseq[i] = '4';
				break;
			case 46://'.'
				cseq[i] = '0';
				break;
			default:
				cseq[i] = '9';
				break;
			}
		}

		String strSeq = new String(cseq);

		int number = strToInt(strSeq);
		return number;
	}

	public int codeFilter(String s) {
		List<String> filters = Arrays.asList(s.split(";"));
		char[] cfilt = new char[filters.size()];
		for (int i = 0; i < filters.size(); i++) {
			String filt = filters.get(i);
			switch (filt.toUpperCase()) {
			case "PASS":
				cfilt[i] = '9';
				break;
			case "AMB":
				cfilt[i] = '1';
				break;
			case "LOWCOV":
				cfilt[i] = '5';
				break;
			case "DEL":
				cfilt[i] = '3';
				break;
			default:
				cfilt[i] = '0';
				break;
			}
		}
		return Integer.parseInt(new String(cfilt));
	}

	public static int strToInt(String str) {
		int i = 0;
		int num = 0;
		boolean isNeg = false;

		// Check for negative sign; if it's there, set the isNeg flag
		if (str.charAt(0) == '-') {
			isNeg = true;
			i = 1;
		}

		// Process each character of the string;
		while (i < str.length()) {
			num *= 10;
			num += str.charAt(i++) - '0'; // Minus the ASCII code of '0' to get
			// the value of the charAt(i++).
		}

		if (isNeg)
			num = -num;
		return num;
	}
/*
	public void printVCFmatrix() {
		System.out.println("Destination path of outputMatrixFile: "+outputMatrixFile);
		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		try {
			myConsole = new PrintStream(new File(outputMatrixFile));
			System.setOut(myConsole);
			for (int i = 0; i < vcfMatrix.size(); i++) {
				System.out.println(vcfMatrix.get(i));
			}
			myConsole.close();
		}  catch (Exception e) {
			System.err.println("Error trying to write outputMatrixFile ");
			e.printStackTrace();
		}  finally {
			if (myConsole != null) {	        	
				myConsole.close();
				System.setOut(stdout);  
			}
		}
		System.setOut(stdout);  
	}
*/
	public void printHumanLines() throws FileNotFoundException {

		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		
		try {
			myConsole = new PrintStream(new File(outputCleanFile));
			System.setOut(myConsole);
			for (int i = 0; i < humanLines.size(); i++) {
				System.out.println(humanLines.get(i));
			}
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
}
