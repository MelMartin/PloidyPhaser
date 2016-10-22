
public class VariationData {
	String id="";//pos+ either(reference or alternative allele)+ '-' before or after indicates an inseert or a deletion
	String name="";
	int pos=0;
	String ref="";
	String alt=".";
	String sample="0/0";
	boolean isInsert=false;
	boolean isDeletion=false;
	int matrixIndex;

	public VariationData(){
		
	}
	public VariationData(String n,int p, String r,String a,String s ){
		name=n;
		pos=p;
		ref=r;
		alt=a;
		sample=s;
		id=Integer.toString(pos)+ref;//set regular id with ref
		if (ref.length()>alt.length())isDeletion=true;
		if (alt.length()>ref.length()){
			isInsert=true;
			id+="-";//set id with insert identifier '-'
		}
	}
	
	public void printOut(){
		System.out.println(this.outString());
	}
	
	public String outString(){
		return /*name+"\t"+*/pos+"\t"+ref+"\t"+alt+"\t"+sample;		
	}
	
	public VariationData makeAlterativeVarDataFromRef(VariationData vd) {
		VariationData nvd=new VariationData(vd.name,vd.pos,vd.ref,vd.alt,vd.sample);
		nvd.id=Integer.toString(pos)+alt;//set regular id with alt
		if (ref.length()>alt.length()){
			isDeletion=true;
			nvd.id+="-";//set id with insert identifier '-'
		}
		return nvd;
	}
	
	public void setMatrixIndex(int i){
		matrixIndex=i;
	}
}
