
public class VariationData {
	String name="";
	int pos=0;
	String ref="";
	String alt=".";
	String sample="0/0";

	public VariationData(){
		
	}
	public VariationData(String n,int p, String r,String a,String s ){
		name=n;
		pos=p;
		ref=r;
		alt=a;
		sample=s;
	}
	
	public void printOut(){
		System.out.println(this.outString());
	}
	
	public String outString(){
		return name+"\t"+pos+"\t"+ref+"\t"+alt+"\t"+sample;		
	}
}
