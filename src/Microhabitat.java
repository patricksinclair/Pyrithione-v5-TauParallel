import org.apache.commons.math3.distribution.LogNormalDistribution;
import java.util.ArrayList;
import java.util.Random;

public class Microhabitat {

    private double mu = Math.log(7.92016113), sigma = 0.10018864;
    private LogNormalDistribution MIC_distribution = new LogNormalDistribution(mu, sigma);

    Random rand = new Random();

    private int K; //karrying kapacity
    private double c; //conc of antimicrobial
    private double b; //migraiton rate
    private ArrayList<Double> population; //list of MICs of bacteria

    private boolean surface = false, biofilm_region = false;
    private double threshold_stickiness; //fraction of occupation required for biofilm classification


    public Microhabitat(int K, double c){
        this.K = K;
        this.c = c;
        this.b = 0.1;
        this.threshold_stickiness = 0.15;
        this.population = new ArrayList<>(K);
    }


    public int getN(){return population.size();}
    public double getC(){return c;}
    public boolean getBiofilm_region(){return this.biofilm_region;}
    public boolean getSurface(){return this.surface;}
    public double getThreshold_stickiness(){return threshold_stickiness;}
    public ArrayList<Double> getPopulation(){return population;}

    public void setSurface(boolean surface){this.surface = surface;}
    public void setBiofilm_region(boolean biofilm_region){this.biofilm_region = biofilm_region;}

    public double fractionFull(){
        //double frac = getN()/(double)K;
        //System.out.println("fraction full: N = "+getN()+"\tK = "+K+"\tfrac = "+frac);
        return getN()/(double)K;
    }

    public double stickiness(){
        //returns parameter between 0 and 1.  Closer to 1 = harder for bacteria to migrate out.
        //TODO change this if we want varying stickiness thresholds etc
        //todo should surface stickiness be 1?
        double sticky_alpha = 0.001630935;
        return (getN() < threshold_stickiness*K) ? 0. : Math.exp(sticky_alpha*(getN() - threshold_stickiness*K))-1.;
    }


    public boolean atCapacity(){
        return fractionFull() > threshold_stickiness;
    }


    public double migrate_rate(){
        return b*(1-stickiness());
    }

    public double beta(int index){
        return population.get(index);
    }

    public double phi_c(int index){
        double cB = c/beta(index);
        return 1. - 6.*cB*cB/(5. + cB*cB);
    }

    public double replicationOrDeathRate(int index){
        return (phi_c(index) > 0.) ? phi_c(index)*(1. - getN()/(double)K) : phi_c(index);
    }


    public double getAvgGenotype(){
        if(getN()==0) return 0.;
        else{
            double sum = 0.;
            for(Double geno : population){
                sum += geno;
            }
            return sum/(double)getN();
        }
    }

    public double getStDevOfGenotype(){
        if(getN() <= 1) return 0.;
        else{
            double mean = getAvgGenotype();
            double sumSq = 0.;

            for(Double geno : population){
                sumSq += (geno-mean)*(geno-mean);
            }
            return Math.sqrt(sumSq/(getN()-1.));
        }
    }


    public void addARandomBacterium_x_N(int n_bacteria){
        for(int i = 0; i < n_bacteria; i++){
            population.add(MIC_distribution.sample());
        }
    }

    public void replicateABacterium_x_N(int index, int nReps){
        for(int i = 0; i < nReps; i++){
            population.add(population.get(index));
        }
    }

    public void addABacterium(double MIC){population.add(MIC);}

    public void removeABacterium(int index){population.remove(index);}



}
