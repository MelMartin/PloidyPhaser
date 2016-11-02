import java.util.ArrayList;
import java.util.List;



public class VariationsConnectivityMatrix {
	ReadList[][]varExpMat;
	
	
	public VariationsConnectivityMatrix (int d) {
		varExpMat=new ReadList[d][d];
		for (int r=0;r<d;r++){
			for (int c=0;c<d;c++){
				varExpMat[r][c]=new ReadList();
			}
		}
	}
	


	public void addConection(int row, int col, Integer readIndex) {
		System.out.println (" VariationsConnectivityMatrix(): ENTERED");
		if(varExpMat[row][col].isNull()){
			 

			varExpMat[row][col].initialize();
		}
		varExpMat[row][col].add(readIndex);
		
	}
	
	
	
	
    private class ReadList{ 
    	
    	List <Integer> ReadsIndexes;
    	
    	private ReadList (){
    		ReadsIndexes=null;
    	}
    	
    	private boolean isNull(){    		
    		System.out.println (" isNull(): PASSED");
    		return ReadsIndexes==null;
    	}
    	
    	private void initialize(){
    		ReadsIndexes=new ArrayList<Integer>();
    	}
    	private boolean add(Integer readIndex){
    		return ReadsIndexes.add(readIndex);
    	}
    	
    }
	
}
