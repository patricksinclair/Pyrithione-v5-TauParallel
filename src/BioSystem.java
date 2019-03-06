import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

public class BioSystem {

    Random rand = new Random();

    private int L, K; //length, karrying kapacity
    private double alpha, c_max; //steepness of antimicrobial gradient, max concn
    private Microhabitat[] microhabitats;
    private double timeElapsed;
    private double tau; //timestep used in tau-leaping
    private double immigration_rate = 100.;

    public BioSystem(int L, int K, double alpha, double c_max, double tau){
        this.L = L;
        this.K = K;
        this.alpha = alpha;
        this.c_max = c_max;
        this.tau = tau;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;

        for(int i = 0; i < L; i++) {
            microhabitats[i] = new Microhabitat(K, BioSystem.getCValWithOffset(i, c_max, alpha, L));
        }

        microhabitats[L - 1].setSurface(true);
        microhabitats[L - 1].setBiofilm_region(true);
        microhabitats[L - 1].addARandomBacterium_x_N(25);
    }


    public double getTimeElapsed(){
        return timeElapsed;
    }

    public int getPop_i(int index){
        return microhabitats[index].getN();
    }

    public int getN_i(int index){
        return microhabitats[index].getN();
    }

    public int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }


    public int getBiofilmEdge(){
        int edgeIndex = 0;
        for(int i = L - 1; i > 0; i--) {
            if(microhabitats[i].getBiofilm_region() && !microhabitats[i - 1].getBiofilm_region()) {
                edgeIndex = i;
                break;
            }
        }
        return edgeIndex;
    }


    public int getBiofilmSize(){
        return L - getBiofilmEdge();
    }


    public double[] getPopulationDistribution(){
        double[] popSizes = new double[L];
        for(int i = 0; i < L; i++) {
            popSizes[i] = microhabitats[i].getN();
        }
        return popSizes;
    }


    public double[] getAvgGenotypeDistribution(){
        double[] avgGenos = new double[L];
        for(int i = 0; i < L; i++) {
            avgGenos[i] = microhabitats[i].getAvgGenotype();
        }
        return avgGenos;
    }


    public double[] getStDevOfGenotypeDistribution(){
        double[] genoStDevs = new double[L];
        for(int i = 0; i < L; i++) {
            genoStDevs[i] = microhabitats[i].getStDevOfGenotype();
        }
        return genoStDevs;
    }


    public void immigrate(int mh_index, int n_immigrants){
        microhabitats[mh_index].addARandomBacterium_x_N(n_immigrants);
    }


    public void migrate(int mh_index, int n_migrations, int originalPopSize){

        int biof_edge = getBiofilmEdge();

        if(originalPopSize > 0) {
            for(int i = 0; i < n_migrations; i++) {

                int rand_bac_index = rand.nextInt(originalPopSize);
                double rand_bac = microhabitats[mh_index].getPopulation().get(rand_bac_index);

                if(mh_index == L - 1) {
                    microhabitats[mh_index].removeABacterium(rand_bac_index);
                    microhabitats[mh_index - 1].addABacterium(rand_bac);
                } else if(mh_index == biof_edge) {
                    //TODO if the bacteria is at the edge of the biofilm, there's a chance it detaches, depending on the stickiness
                    microhabitats[mh_index].removeABacterium(rand_bac_index);
                    microhabitats[mh_index + 1].addABacterium(rand_bac);
                } else {
                    if(rand.nextBoolean()) {
                        microhabitats[mh_index].removeABacterium(rand_bac_index);
                        microhabitats[mh_index + 1].addABacterium(rand_bac);
                    } else {
                        microhabitats[mh_index].removeABacterium(rand_bac_index);
                        microhabitats[mh_index - 1].addABacterium(rand_bac);
                    }
                }
            }
        }
    }


    public void updateBiofilmSize(){
        //todo change the constraint on i if issues arise
        int i_lim = Math.max(getBiofilmEdge() - 1, 1);
        for(int i = L - 1; i >= i_lim; i--) {
            if(microhabitats[i].atCapacity()) {
                microhabitats[i].setBiofilm_region(true);

                if(!microhabitats[i - 1].atCapacity()) {
                    microhabitats[i - 1].setBiofilm_region(true);
                }
            }
        }
    }


    public void performAction(){

        Microhabitat[] updated_microhabs = microhabitats.clone();

        double tau_step = tau;
        int biofilm_size = getBiofilmSize();
        int[][] replicationAllocations;
        int[][] deathAllocations;
        int[] migrationAllocations;
        int n_immigrants;

        whileloop:
        while(true) {

            replicationAllocations = new int[biofilm_size][];
            deathAllocations = new int[biofilm_size][];
            migrationAllocations = new int[biofilm_size];

            int mh_counter = 0; //this loop determines the number of events for each bacteria/microhab
            for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++) {

                int mh_pop = microhabitats[mh_index].getN();
                int[] n_replications = new int[mh_pop];
                int[] n_deaths = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++) {

                    double g_or_d_rate = microhabitats[mh_index].replicationOrDeathRate(bac_index);

                    if(g_or_d_rate == 0.) {
                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = 0;

                    } else if(g_or_d_rate > 0) {
                        n_replications[bac_index] = new PoissonDistribution(g_or_d_rate*tau_step).sample();
                        n_deaths[bac_index] = 0;

                    } else {
                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = new PoissonDistribution(Math.abs(g_or_d_rate)*tau_step).sample();
                        if(n_deaths[bac_index] > 1) {
                            tau_step /= 2.;
                            //System.out.println("DOUBLE DEATH HANDLED");
                            continue whileloop;
                        }
                    }
                }

                double migration_rate = microhabitats[mh_index].migrate_rate();
                int n_migrations = (migration_rate > 0) ? new PoissonDistribution(migration_rate*tau_step).sample() : 0;

                replicationAllocations[mh_counter] = n_replications;
                deathAllocations[mh_counter] = n_deaths;
                migrationAllocations[mh_counter] = n_migrations;

                mh_counter++;
            }
            n_immigrants = new PoissonDistribution(immigration_rate*tau_step).sample();
            break whileloop;
        }


        int mh_counter_growth = 0; //this loop carries out the calculated actions for replication and death
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++) {

            int originalPopSize = microhabitats[mh_index].getN();

            for(int bac_index = originalPopSize - 1; bac_index >= 0; bac_index--) {
                updated_microhabs[mh_index].replicateABacterium_x_N(bac_index, replicationAllocations[mh_counter_growth][bac_index]);
                if(deathAllocations[mh_counter_growth][bac_index] > 0)
                    updated_microhabs[mh_index].removeABacterium(bac_index);
            }
            mh_counter_growth++;
        }


        int mh_counter_migration = 0;
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++) {

            int originalPopSize = microhabitats[mh_index].getN();
            migrate(mh_index, migrationAllocations[mh_counter_migration], originalPopSize);
            mh_counter_migration++;
        }

        microhabitats = updated_microhabs;
        immigrate(getBiofilmEdge(), n_immigrants);
        //System.out.println("n_immigrants: "+n_immigrants);
        updateBiofilmSize();
        timeElapsed += tau_step;
    }


    public static void tester(){

        double duration = 200;
        int L = 500;
        int K = 500;
        double c_max = 10.;
        double alpha = 0.01;
        double tau = 0.01;
        int r = 0;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

        while(bs.getTimeElapsed() <= duration) {

            String output = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tbf_edge pop %d \tbf_edge fracfull: %.3f",
                    r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bs.getBiofilmEdge()), bs.microhabitats[bs.getBiofilmEdge()].fractionFull());

            String output2 = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d",
                    r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge());

            System.out.println(output);

            bs.performAction();
        }

    }


    public static void getPopDistbInfo(){
        //method to get info on population distbs
        //get popsize over time
        //pop distb over time
        //biofilm edge over time
        //avg genotype distb over time

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nReps = 10, nMeasurements = 20;
        double duration = 100., interval = duration/nMeasurements;

        String popSizeFilename = "pyrithione-testing-pop_size-t=" + String.valueOf(duration);
        String popDistbFilename = "pyrithione-testing-pop_distb-t=" + String.valueOf(duration);
        String biofilmEdgeFilename = "pyrithione-testing-biofilm_edge-t=" + String.valueOf(duration);
        String avgGenotypeDistbFilename = "pyrithione-testing-avgGenoDistb-t=" + String.valueOf(duration);
        String genoStDevDistbFilename = "pyrithione-testing-genoStDevDistb-t=" + String.valueOf(duration);
        String counterDistbsFilename = "pyrithione-testing-counterDistb-t=" + String.valueOf(duration);

        String[] counterHeaders = {"immigration", "migrationIn", "migrationOut", "replication", "death"};
        int nCounters = counterHeaders.length;

        double[][] allPopSizes = new double[nReps][];
        double[][][] allPopDistbs = new double[nReps][][];
        double[][] allBiofilmEdges = new double[nReps][];
        double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        double[][][] allGenoStDevs = new double[nReps][][];
        double[][] allCounters = new double[nReps][];


        for(int r = 0; r < nReps; r++) {

            BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);
            //double[] cVals = bs.getCVals();
            //System.out.println(Arrays.toString(cVals));
            //System.out.println(cVals[476]);

            boolean alreadyRecorded = false;
            int timerCounter = 0;

            double[] popSizes = new double[nMeasurements + 1];
            double[][] popDistbs = new double[nMeasurements + 1][];
            double[] biofilmEdges = new double[nMeasurements + 1];
            double[][] avgGenotypeDistbs = new double[nMeasurements + 1][];
            double[][] genoStDevs = new double[nMeasurements + 1][];


            while(bs.timeElapsed <= duration + 0.02*interval) {

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.1*interval) && !alreadyRecorded) {

                    int bf_edge_i = bs.getBiofilmEdge();

                    String output = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tbf_edge pop %d \tbf_edge fracfull: %.3f",
                            r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bf_edge_i), bs.microhabitats[bf_edge_i].fractionFull());

                    String output2 = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d, bf_edge pop: %d",
                            r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bf_edge_i));

                    System.out.println(output2);

                    popSizes[timerCounter] = bs.getTotalN();
                    popDistbs[timerCounter] = bs.getPopulationDistribution();
                    biofilmEdges[timerCounter] = bs.getBiofilmEdge();
                    avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                    genoStDevs[timerCounter] = bs.getStDevOfGenotypeDistribution();

                    alreadyRecorded = true;
                    timerCounter++;
                }
                if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

                bs.performAction();
            }

            allPopSizes[r] = popSizes;
            allPopDistbs[r] = popDistbs;
            allBiofilmEdges[r] = biofilmEdges;
            allAvgGenotypeDistbs[r] = avgGenotypeDistbs;
            allGenoStDevs[r] = genoStDevs;
        }

        double[] processedPopSizes = Toolbox.averagedResults(allPopSizes);
        double[][] processedPopDistbs = Toolbox.averagedResults(allPopDistbs);
        double[] processedBiofilmEdges = Toolbox.averagedResults(allBiofilmEdges);
        double[][] processedAvgGenotypeDistbs = Toolbox.averagedResults(allAvgGenotypeDistbs);
        double[][] processedGenoStDevs = Toolbox.averagedResults(allGenoStDevs);

        Toolbox.writeAveragedArrayToFile(popSizeFilename, processedPopSizes);
        Toolbox.writeAveragedDistbsToFile(popDistbFilename, processedPopDistbs);
        Toolbox.writeAveragedArrayToFile(biofilmEdgeFilename, processedBiofilmEdges);
        Toolbox.writeAveragedDistbsToFile(avgGenotypeDistbFilename, processedAvgGenotypeDistbs);
        Toolbox.writeAveragedDistbsToFile(genoStDevDistbFilename, processedGenoStDevs);

        System.out.println("results written to file");
    }



    public static double[][] getAvgGenoDistbs(int i){

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nMeasurements = 20;
        double duration = 100., interval = duration/nMeasurements;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

        boolean alreadyRecorded = false;
        int timerCounter = 0;

        double[][] avgGenotypeDistbs = new double[nMeasurements + 1][];

        while(bs.timeElapsed <= duration + 0.02*interval) {

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.01*interval) && !alreadyRecorded) {

                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed());
                avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                alreadyRecorded = true;
                timerCounter++;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        return avgGenotypeDistbs;
    }


    public static DataBox getAllData(int i){

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nMeasurements = 40;
        double duration = 200., interval = duration/nMeasurements;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

        boolean alreadyRecorded = false;
        int timerCounter = 0;

        double[] popSizes = new double[nMeasurements+1];
        double[][] popDistbs = new double[nMeasurements+1][];
        double[] biofilmEdges = new double[nMeasurements+1];
        double[][] avgGenotypeDistbs = new double[nMeasurements+1][];
        double[][] genoStDevs = new double[nMeasurements+1][];

        while(bs.timeElapsed <= duration+0.02*interval){

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getBiofilmSize()*K;
                int total_N = bs.getTotalN();

                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+total_N+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge());

                popSizes[timerCounter] = total_N;
                popDistbs[timerCounter] = bs.getPopulationDistribution();
                biofilmEdges[timerCounter] = bs.getBiofilmEdge();
                avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                genoStDevs[timerCounter] = bs.getStDevOfGenotypeDistribution();

                alreadyRecorded = true;
                timerCounter++;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        return new DataBox(popSizes, popDistbs, biofilmEdges, avgGenotypeDistbs, genoStDevs);
    }


    public static void getInfoInParallel(){

        long startTime = System.currentTimeMillis();

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nReps = 16;
        double duration = 200.;

        double[][] allPopSizes = new double[nReps][];
        double[][][] allPopDistbs = new double[nReps][][];
        double[][] allBiofilmEdges = new double[nReps][];
        double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        double[][][] allGenoStDevs = new double[nReps][][];

        String popSizeFilename = "pyrithione-testing-pop_size-t="+String.valueOf(duration)+"-parallel";
        String popDistbFilename = "pyrithione-testing-pop_distb-t="+String.valueOf(duration)+"-parallel";
        String biofilmEdgeFilename = "pyrithione-testing-biofilm_edge-t="+String.valueOf(duration)+"-parallel";
        String avgGenotypeDistbFilename = "pyrithione-testing-avgGenoDistb-t="+String.valueOf(duration)+"-parallel";
        String genoStDevDistbFilename = "pyrithione-testing-genoStDevDistb-t="+String.valueOf(duration)+"-parallel";

        //double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        DataBox[] dataBoxes = new DataBox[nReps];

        IntStream.range(0, nReps).parallel().forEach(i -> dataBoxes[i] = BioSystem.getAllData(i));

        for(int j = 0; j < dataBoxes.length; j++){
            allPopSizes[j] = dataBoxes[j].getPopSizes();
            allPopDistbs[j] = dataBoxes[j].getPopDistbs();
            allBiofilmEdges[j] = dataBoxes[j].getBiofilmEdges();
            allAvgGenotypeDistbs[j] = dataBoxes[j].getAvgGenotypeDistbs();
            allGenoStDevs[j] = dataBoxes[j].getGenoStDevs();
        }

        double[] processedPopSizes = Toolbox.averagedResults(allPopSizes);
        double[][] processedPopDistbs = Toolbox.averagedResults(allPopDistbs);
        double[] processedBiofilmEdges = Toolbox.averagedResults(allBiofilmEdges);
        double[][] processedAvgGenotypeDistbs = Toolbox.averagedResults(allAvgGenotypeDistbs);
        double[][] processedGenoStDevs = Toolbox.averagedResults(allGenoStDevs);

        Toolbox.writeAveragedArrayToFile(popSizeFilename, processedPopSizes);
        Toolbox.writeAveragedDistbsToFile(popDistbFilename, processedPopDistbs);
        Toolbox.writeAveragedArrayToFile(biofilmEdgeFilename, processedBiofilmEdges);
        Toolbox.writeAveragedDistbsToFile(avgGenotypeDistbFilename, processedAvgGenotypeDistbs);
        Toolbox.writeAveragedDistbsToFile(genoStDevDistbFilename, processedGenoStDevs);

        long finishTime = System.currentTimeMillis();

        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);

        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }







    public static double getCValWithOffset(int index, double maxC, double alpha, int L){
        //this calculates i* for the gradient profile offset, moves so the final concn is maxC, and it decreases with 1/e
        //or something like that
        //then calculates the corresponding concentration in that microhabitat

        double offset =  (L-1.) - Math.log(maxC+1.)/alpha;
        return (index >= offset) ? Math.exp(alpha*(index - offset)) - 1. : 0.;
    }
}
