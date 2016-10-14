
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
	List<String> cleanLines = new ArrayList<String>();
	String outputCleanFile="" ;
	boolean preventDoubleLines=false;
	File vcfFile;
	HashMap<VariationData,Integer> variations=new HashMap<VariationData,Integer>();//keeps track of vcf positions<key> 
	//that contains a 'ref' different than an 'alt', and their Index order <value>,in order to build the matrix,
	int[][]varMat;
	//----------------------- AM I GOING TO NEED THESE????--------------------------
	List<Integer> refInsertsPos=new ArrayList<Integer>();//keeps track of inserts pos into the reference allele (each time 'alt' is longer than 'ref')
	List<Integer> nbInsertsPos=new ArrayList<Integer>();//keeps track of nb of inserts into the reference allele (each time 'alt' is longer than 'ref')
	//-----------------------ERASE IF NOT USED -------------------------------------

	
	public VCFManager(File vf) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		//solve the paths
		vcfFile=vf;
		//vcfVariantExtractor(vcfFile.getAbsolutePath());
		vcfVariantShortExtractor(vcfFile.getAbsolutePath());//extracts only the repeated variants in a short format: chromName, pos, ref, alt
		//vcfExtractor(vcfFile.getAbsolutePath());	
		fillVariantMatrix();
		
	}
	private void fillVariantMatrix() {
		varMat=new int[variations.size()][variations.size()];
		SAMFileReader inputSam;
		for (int bi=0;bi<PloidyPhaser.bamPaths.size();bi++){
			inputSam = new SAMFileReader(PloidyPhaser.bamPaths.get(bi));
			
			
			
			
			
			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter = inputSam.iterator();
			//PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//"+Ploest.projectName+"SamParsed.txt", "UTF-8");
			//String line = "";
			int ct = 0;
			//System.out.println("Analyzing ");
			String refName = "";// for debugging a bad line in the .sam file
			int alStart;
			int alEnd;
			

			while (iter.hasNext()) {//iterates the sam file

				SAMRecord rec = iter.next();
				refName = rec.getReferenceName();
				alStart = rec.getAlignmentStart();
				alEnd = rec.getReadLength()+alStart;
				if (alStart<10000 && ct<10 ){
					System.out.println("ct:"+ct+" read:"+rec.getReadName()+" refName:"+refName+" st:"+alStart+" end:"+alEnd+" length:"+rec.getReadLength());
					ct++;
				}
				//line = (i + " " + rec.getReadName() + " " + refName + " " + alStart + " " + rec.getAlignmentEnd() + " "
				//			+ rec.getReadLength() + " " + rec.getMappingQuality() + " ;");
				//writer.println(line);

				
				

			}
			iter.close();
			inputSam.close();
			
		}
		
		
	}
	//extracts a new vcf file with no double lines, and the posibility to modify the sample field (ie:0/0/0)
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
		String prevAlt=".";//alternative allele when existing in a double line
		int qual=0;
		List<String> infoFields=new ArrayList();

		String format="GT";
		String oldSample="0/0";
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
				//oldSample=sc.next();//if we want to keep the format that HPoP writer says then turn on this two comments
				//sample="";//(this one too) plus change inside loop sample for oldSample and turn on the line after the  if
				sample=sc.next();
				
				if( ploidy>2){					
					for (int p=2;p<ploidy;p++){
						sample+="/0";//adds 0/ to the bi-allelic sample field from the left to fill the ploidies >2
					}
				}
				//sample+=oldSample;//(turn this ne on too)
				//if (ct<20)System.out.println(sample);
				
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
				if(preventDoubleLines){
					if(pos!=prevPos && prevPos!=0){//prevents double lines for same position
						cleanLines.add(previousHumanLine);	
					}/*else{
					if(!prevAlt.equals(".")){
						System.out.println("ct:"+ct+" "+previousHumanLine);
						System.out.println("   "+ct+" "+currentHumanLine);
						}
					}
				*/
				}
					else{
						if(!previousHumanLine.equals( "//>>"))cleanLines.add(previousHumanLine);
					}
				
				prevPos=pos;
				prevAlt=alt;
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
			if(preventDoubleLines){
				printCleanLines("Clean");
			}else printCleanLines("");
			

			if (sc != null)
				sc.close();
			
		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentHumanLine);
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

	
//creates a new vcf file with only the variants on it (SNPs, structural variations and double lines) 	
public void vcfVariantExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

		
		
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
		Boolean isVariation=false;
		int ploidy=0;
		String ref = "";// reference allele
		String alt=".";// alternative allele
		String prevAlt=".";//alternative allele when existing in a double line
		int qual=0;
		List<String> infoFields=new ArrayList();

		String format="GT";
		String oldSample="0/0";
		String sample="";
		String nextFilter;
		String currentHumanLine = "*";
		String previousHumanLine = "//>>";
		String lastWritenLine="";

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
				//oldSample=sc.next();//if we want to keep the format that HPoP writer says then turn on this two comments
				//sample="";//(this one too) plus change inside loop sample for oldSample and turn on the line after the  if
				sample=sc.next();
				
				if( ploidy>2){					
					for (int p=2;p<ploidy;p++){
						sample+="/0";//adds 0/ to the bi-allelic sample field from the left to fill the ploidies >2
					}
				}
				//sample+=oldSample;//(turn this one on too)
				
				
				if (!infoFields.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )

					svLength = Integer.parseInt(infoFields.get(1).substring(6, infoFields.get(1).length())) ;
					end=Integer.parseInt(infoFields.get(2).substring(4, infoFields.get(2).length())) ;
					
					System.out.println("pos:"+pos+" SV infoFields:"+infoFields+ " LENGTH:" + svLength+" FORMAT:"+format);
					currentHumanLine = chrom+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" +qual+"\t"+ nextFilter + "\tSVTYPE="
							+ infoFields.get(0) + ";SVLEN=" + svLength+";END="+end+"\t"+format+"\t"+sample;
					isVariation=true;
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
				
				if(pos==prevPos && !lastWritenLine.equals(previousHumanLine)){
						cleanLines.add(previousHumanLine);
						lastWritenLine=previousHumanLine;
					}
				if((isVariation || !alt.equals("."))&& !lastWritenLine.equals(currentHumanLine)){
					cleanLines.add(currentHumanLine);
					isVariation=false;
					lastWritenLine=currentHumanLine;
				}
				prevPos=pos;
				prevAlt=alt;
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
		
			printCleanLines("Repeats");

			if (sc != null)
				sc.close();
			
		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentHumanLine);
		}

	}


public void vcfVariantShortExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

	
	
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
	Boolean isVariation=false;
	int ploidy=0;
	String ref = "";// reference allele
	String alt=".";// alternative allele
	String prevAlt=".";//alternative allele when existing in a double line
	int qual=0;
	List<String> infoFields=new ArrayList();

	String format="GT";
	String oldSample="0/0";
	String sample="";
	String nextFilter;
	String currentHumanLine = "*";
	String previousHumanLine = "//>>";
	String lastWritenLine="";

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
		while (sc.hasNextLine()  && pos<10010 /* */) {//test with only the first 10.000 positions
			
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
			//----------------------------------------------SAMPLE FIELD: Two possibilities; The current one works, but HPoP writer says that the other is needed for its algorythm to work
			//oldSample=sc.next();//if we want to keep the format that HPoP writer says,('0/'+'0/1'), then turn on this + two indicated comments in this same block
			//sample="";//(this one too) plus change inside loop :'sample' for 'oldSample'
			sample=sc.next();
			
			if( ploidy>2){					
				for (int p=2;p<ploidy;p++){
					sample+="/0";//adds '/0' to the bi-allelic sample field from the right('0/1'+'/0') to fill the ploidies >2
								//if we want to keep the format that HPoP writer says,  change :'sample' for 'oldSample' and '/0' for '0/'
				}
			}
			//sample+=oldSample;//(turn this one on too)if we want to keep the format that HPoP writer says
			//----------------------------------------------
			
			//buid the line that is going to be written
			if (!infoFields.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )

				svLength = Integer.parseInt(infoFields.get(1).substring(6, infoFields.get(1).length())) ;
				end=Integer.parseInt(infoFields.get(2).substring(4, infoFields.get(2).length())) ;
				
				System.out.println("pos:"+pos+" SV infoFields:"+infoFields+ " LENGTH:" + svLength+" FORMAT:"+format);
				currentHumanLine = chrom+"\t" + pos + "\t"+ ref + "\t" + alt +"\t"+sample;
				isVariation=true;
			} else {// ...regular vcf line
				svLength = 0;
				//System.out.print("ct:"+ct+"\t" + pos + "\t.\t"+ ref + "\t" + alt + "\t" );
				currentHumanLine = chrom+"\t" + pos + "\t"+ ref + "\t" + alt + "\t"+sample ;		
			}
			
			//write out the builded lines
			/*
			if(pos==prevPos && !lastWritenLine.equals(previousHumanLine) ){//if is a repeated position (but different entry)
					cleanLines.add("= "+previousHumanLine);
					lastWritenLine=previousHumanLine;
				}
			*/
			if((isVariation || !alt.equals("."))/*&& !lastWritenLine.equals(currentHumanLine)*/){
				cleanLines.add("* "+currentHumanLine);
				isVariation=false;
				//lastWritenLine=currentHumanLine;
				variations.put(new VariationData(pos, ref, alt),(variations.size()+1));
System.out.println("pos:"+pos+"/"+ref+"/"+alt+" ind:"+variations.size());
				if (alt.length()>ref.length()){
					refInsertsPos.add(pos);
					nbInsertsPos.add(alt.length()-ref.length());
				}
			}
			prevPos=pos;
			prevAlt=alt;
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
	
		printCleanLines("Short");

		if (sc != null)
			sc.close();
		
	} catch (Exception e) {
		System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentHumanLine);
	}

}
	
}
