
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;


public class VCFParser {
	
	String header="";//vcf Header
	
	//List<VariationData> variations = new ArrayList<VariationData>();
	String outputCleanFile="" ;
	File vcfFile;
	VariationsManager varMan=new VariationsManager(this);
	
	
	public VCFParser(File vf) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		
		vcfFile=vf; 
		vcfVariantShortExtractor(vcfFile.getAbsolutePath());//extracts only the repeated variants in a short format: chromName, pos, ref, alt
		varMan.constructVariantMatrix();
		varMan.fillVariantMatrix();
		varMan.colourMatrix(2);//activate vcm.printMatrix(20); in varMan.fillVariantMatrix()

	}

	

	//extracts the vcf lines that contains a variation with a ref different than an alternative allele (both !='.')
	public void vcfVariantShortExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

		String line="";
		Scanner sc = new Scanner(new File(inputFile));

		int ct = 0;
		String chrom="";
		int pos = 0;
		int svLength;// length of the sv
		//int end;
		Boolean isPloidySolved=false;
		Boolean isVariation=false;
		int ploidy=0;
		String ref = "";// reference allele
		String alt=".";// alternative allele
		List<String> infoFields=new ArrayList<String>();

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
			//solve sample field (general for all lines)
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
					//end=Integer.parseInt(infoFields.get(2).substring(4, infoFields.get(2).length())) ;
					System.out.println("pos:"+pos+" SV infoFields:"+infoFields+ " LENGTH:" + svLength+" FORMAT:"+format);
					currentVariationData = new VariationData(chrom, pos ,ref,alt,sample);
					isVariation=true;
					//svLength = 0;
				} else {// ...regular vcf line
					currentVariationData =new VariationData(chrom, pos ,ref,alt,sample);		
				}
				

				//write out the builded lines (these are the ones that interest us)
				if((isVariation || !alt.equals("."))){
					varMan.vcfVarLines.add(currentVariationData);//adds each vcf line with ref!=alt (
					isVariation=false;
					varMan.registerPossibleVariationExpressions(currentVariationData);//register separately each of the possible expressions per variation line ('ref' and 'alt')
				}

				if (sc.hasNextLine() ) { // 'if' to avoid error at end of file
					currentVariationData =  new VariationData();
					chrom=sc.next();// contig name 
				}else {
					varMan.vcfVarLines.add(currentVariationData);//adds last line
				}
				ct++;		
			}

			varMan.printOutVariations("Short");

			if (sc != null)
				sc.close();
			
			//CreateVariationVertex();//creates the vertex that are going to be used to define the phasing Matrix 'varMat'

		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentVariationData.outString());
		}
	}


	
}
