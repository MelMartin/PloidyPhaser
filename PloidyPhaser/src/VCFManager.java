
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;


public class VCFManager {
	List<String> cleanLines = new ArrayList<String>();
	String outputCleanFile="" ;


	
	public VCFManager(File vcfFile) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		//solve the paths
		System.out.println("outputCleanFile:"+vcfFile.getParent()+"\\"+vcfFile.getName().substring(0, vcfFile.getName().lastIndexOf('.'))+"Clean.vcf");
		outputCleanFile = vcfFile.getParent()+"\\"+vcfFile.getName().substring(0, vcfFile.getName().lastIndexOf('.'))+"Clean.vcf";
		
		vcfExtractor(vcfFile.getAbsolutePath());	
	}
	
	public void vcfExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

		
		
		String line="";
		String header="";
		Scanner sc = new Scanner(new File(inputFile));
		
		int ct = 0;
		String chrom="";
		int prevPos=0;
		int pos = 0;

		int svLength;// length of the sv
		int end;
		Boolean isPloidySolved=false;
		int ploidy=0;
		String ref = "";// reference allele
		String alt=".";// alternative allele
		int qual=0;
		List<String> infoFields=new ArrayList();

		String format="GT";
		String sample="";
		String nextFilter;
		String currentHumanLine = "*";
		String previousHumanLine = "//>>";

		String next="";
		try {
			// skip the header		
			line=sc.nextLine();
			while (line.substring(0, 1).equals("#")){				
				header+=line ;	
				line=sc.nextLine();				
				if(line.substring(0, 1).equals("#"))header+="\n";
				ct++;
			}

			cleanLines.add(header);
			
			//reset scaner after the header;
			sc = new Scanner(new File(inputFile));
			for(int l=0;l<ct;l++){
				sc.nextLine();
			}
			ct=0;
			chrom=sc.next();
		
			//get the values
			while (sc.hasNextLine() /* && ct < 34000 */) {
				
				pos = Integer.parseInt(sc.next());	// get pos	
				sc.next();// skip id
				ref = sc.next();// get reference allele
				alt = sc.next();// get alternative allele
				next=sc.next();//qual
				
				if(!next.equals(".")){
					//System.out.print("-ct:"+ct+" pos:"+pos);
					qual=Integer.parseInt(next);//  qual
					//System.out.println(" qual:"+qual);
					
				}else {
					//System.out.println("-.ct:"+ct+" pos:"+pos);
					qual=0;
					//System.out.println( " qual == '.' :"+qual);
				}
				//
				nextFilter = sc.next();
				
				infoFields = Arrays.asList(sc.next().split(";"));// split all
				// fields in the info line
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

					svLength = Integer.parseInt(infoFields.get(1).substring(6, infoFields.get(1).length())) ;
					end=Integer.parseInt(infoFields.get(2).substring(4, infoFields.get(2).length())) ;
					
					System.out.println("pos:"+pos+" SV infoFields:"+infoFields+ " LENGTH:" + svLength+" FORMAT:"+format);
					currentHumanLine = chrom+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" +qual+"\t"+ nextFilter + "\tSVTYPE="
							+ infoFields.get(0) + ";SVLEN=" + svLength+";END="+end+"\t"+format+"\t"+sample;

				} else {// ...regular vcf line
					
					svLength = 0;
					//System.out.print("ct:"+ct+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" );
					//System.out.print(qual+"\t"+nextFilter + "\t");
					currentHumanLine = chrom+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" + +qual+"\t"+nextFilter + "\t";
					for (int f=0;f<(infoFields.size()-1);f++){
						//System.out.println("f:"+f+infoFields.get(f));
						currentHumanLine += infoFields.get(f)+";";
					}
					//System.out.println("ct:"+ct+" f:"+(infoFields.size()-1)+infoFields.get(infoFields.size()-1));
					currentHumanLine +=   infoFields.get(infoFields.size()-1)+"\t"+format+ "\t"+sample;
					
				}
				
				if(pos!=prevPos && prevPos!=0){
					cleanLines.add(previousHumanLine);
				}
				prevPos=pos;
				previousHumanLine = currentHumanLine;
				
				//line = sc.nextLine();
				if (sc.hasNextLine() ) { // 'if' to avoid error at end of file
					currentHumanLine = sc.nextLine();
					chrom=sc.next();// contig name 
				}else {
					cleanLines.add(previousHumanLine);//adds last line
					//System.out.println("no next line at pos "+pos);
				}
				ct++;
			
				
			}
		
			printCleanLines();

			if (sc != null)
				sc.close();
			
		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentHumanLine);
		}

	}
	

	public void printCleanLines() throws FileNotFoundException {

		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		
		try {
			myConsole = new PrintStream(new File(outputCleanFile));
			System.setOut(myConsole);
			for (int i = 0; i < cleanLines.size(); i++) {
				System.out.println(cleanLines.get(i));
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

	
}
