/*
package edu.uw.waterlooms.database;

import edu.uw.waterlooms.match.Genome;
import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.match.IsotopeTrail;
import edu.uw.waterlooms.match.Peptide;

import java.io.FileInputStream;
import java.io.IOException;
import java.sql.*;
import java.util.*;

public class JDBCPostgresService {
    private String url;
    private Properties props;
    private Connection conn;

    public JDBCPostgresService(){
        this.url = "jdbc:postgresql://localhost:5432/waterlooms";
        this.props = new Properties();
        this.props.setProperty("user", "waterlooms");
        this.props.setProperty("password", "waterlooms");
        initializeDBConnection();
    }

    private void initializeDBConnection(){
        try {
            this.conn = DriverManager.getConnection(this.url, this.props);
        } catch (Exception $e) {
            System.err.println("Unable to acquire DB Connection.");
            $e.printStackTrace();
        }
    }

    public int insertProtein(Genome $genome) throws SQLException{
        String query = "INSERT INTO protein (protein_header, protein_sequence, protein_description) values(?, ?, ?)";
        String generatedColumns[] = { "id" };
        PreparedStatement statement = this.conn.prepareStatement(query, generatedColumns);
        statement.setString(1, $genome.id);            // protein_header
        statement.setString(2, "");                 // protein_sequence TODO: Finish this
        statement.setString(3, $genome.description);  // protein_description
        statement.execute(); // Commit result to DB

        ResultSet insertedId = statement.getGeneratedKeys();
        if (insertedId.next()){
            return insertedId.getInt(1);
        } else {
            throw new SQLException("Unable to retrieve ID of inserted protein.");
        }
    }

    public int getProteinByHeader(String $proteinHeader) throws SQLException {
        String query = "SELECT id FROM protein where protein_header = ?";
        PreparedStatement statement = this.conn.prepareStatement(query);
        statement.setString(1, $proteinHeader);

        ResultSet rs = statement.executeQuery();
        if (rs.next()){
            return rs.getInt("id");
        } else {
            throw new SQLException("Unable to find protein with header: " + $proteinHeader);
        }
    }

    public void insertIsolationWindowAndTrails(String ms2xicFile) throws IOException, SQLException {
        //ArrayList<IsolationWindow> windows =  new ArrayList();
        String dataFile = ms2xicFile;
        // create new isolation window
        IsolationWindow cur = new IsolationWindow();
        IsotopeTrail curTrail = new IsotopeTrail();
        int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail
        int isolationWindowID = 0;

        FileInputStream inputStream = null;
        Scanner sc = null;
        try {
            inputStream = new FileInputStream(dataFile);
            sc = new Scanner(inputStream, "UTF-8");
            while (sc.hasNextLine()) {
                String line = sc.nextLine();
                if (line.isEmpty()) { continue; }
                if (line.contains("END")){
                    // Sort Trails Here
                    Collections.sort(cur.trails);
                    //windows.add(cur);
                    continue;
                }
                if (line.contains("START")) {

                    cur = new IsolationWindow();
                    line = sc.nextLine();
                    String[] range = line.split("\\s+");
                    cur.mzLow = Double.parseDouble(range[0]);
                    cur.mzHigh = Double.parseDouble(range[1]);

                    // Insert IsolationWindow here
                    String query = "insert into isolationwindow (start_mz, end_mz, start_rt, end_rt)" + " values (?, ?, ?, ?)";
                    String generatedColumns[] = {"id"};
                    PreparedStatement isolationWindowStatement = conn.prepareStatement(query, generatedColumns);
                    isolationWindowStatement.setDouble(1, cur.mzLow);
                    isolationWindowStatement.setDouble(2, cur.mzHigh);
                    isolationWindowStatement.setDouble(3, 0);
                    isolationWindowStatement.setDouble(4, 0);
                    isolationWindowStatement.execute();
                    ResultSet insertedId = isolationWindowStatement.getGeneratedKeys();
                    if (insertedId.next())
                    {
                        isolationWindowID = insertedId.getInt(1);
                    }

                    cur.trails = new ArrayList<>();
                    lineType = 0;
                    continue;
                }
                if (lineType == 0) {
                    String[] mzLine = line.split("\\s+");
                    curTrail.mz = Double.parseDouble(mzLine[0]);
                    curTrail.maxIntensityRt = Double.parseDouble(mzLine[1]);
                    lineType = 1;
                    continue;
                }
                if (lineType == 1) {
                    String[] linelst = line.split("\\s+");
                    curTrail.ParseIntensities(linelst);
                    lineType = 2;
                    continue;
                }
                if (lineType == 2) {
                    String[] linelst = line.split("\\s+");
                    curTrail.ParseRts(linelst);
                    // add to list
                    cur.PutRtMax(curTrail.rts, curTrail.intensities);

                    String query = "insert into trail (isolation_window_id, local_max_mz, local_max_rt, intensities, retention_times) values(?, ?, ?, ?, ?)";
                    PreparedStatement trailStatement = conn.prepareStatement(query);
                    trailStatement.setInt(1, isolationWindowID);
                    trailStatement.setDouble(2, curTrail.mz);
                    trailStatement.setDouble(3, curTrail.maxIntensityRt * 60);

                    Array intensityArray = conn.createArrayOf("float",curTrail.intensities);
                    trailStatement.setArray(4, intensityArray);

                    Array retentionTimeArray = conn.createArrayOf("float",
                            Arrays.stream(curTrail.rts).map(x -> x * 60).toArray()
                    );
                    trailStatement.setArray(5, retentionTimeArray);
                    trailStatement.execute();


                    //cur.trails.add(curTrail);
                    curTrail = new IsotopeTrail();
                    lineType = 0;
                    continue;
                }
                // System.out.println(line);
            }
            // note that Scanner suppresses exceptions
            if (sc.ioException() != null) {
                throw sc.ioException();
            }
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
            if (sc != null) {
                sc.close();
            }
        }
    }


    public void insertCandidatePeptideAndFragmentIons(Peptide $candidatePeptide) {
        try {
            // Perform a lookup to see the protein id
            int protein_id = this.getProteinByHeader($candidatePeptide.id);

            String query = "INSERT INTO candidatepeptide (protein_id, unmodified_sequence, modified_sequence, charge, theoretical_mz) values (?, ?, ?, ?, ?)";
            String generatedColumns[] = {"id"};
            PreparedStatement statement = this.conn.prepareStatement(query, generatedColumns);
            statement.setInt(1, protein_id);                        // protein_header
            statement.setString(2, $candidatePeptide.composition);  // protein_sequence
            statement.setString(3, null);
            statement.setInt(4, $candidatePeptide.charge);                             // charge
            statement.setDouble(5, $candidatePeptide.mz);          // protein_mz
            statement.execute(); // Commit result to DB

            ResultSet insertedId = statement.getGeneratedKeys();
            if (insertedId.next()) {
                int peptide_id = insertedId.getInt(1);

                // Iterate through the Fragment Ions y/b inserting them into the database.
                for (int i = 0; i < $candidatePeptide.y_ions.length; i++) {
                    // Insert the Y-Ion
                    query = "INSERT INTO fragmention (candidate_peptide_id, ion_type, ion_position, theoretical_mz) values (?, ?, ?, ?)";
                    statement = this.conn.prepareStatement(query);
                    statement.setInt(1, peptide_id);
                    statement.setString(2, "y");
                    statement.setInt(3, i + 1);
                    statement.setDouble(4, $candidatePeptide.y_ions[i]);
                    statement.execute();

                    // Insert the B-Ion
                    statement = this.conn.prepareStatement(query);
                    statement.setInt(1, peptide_id);
                    statement.setString(2, "b");
                    statement.setInt(3, i + 1);
                    statement.setDouble(4, $candidatePeptide.b_ions[i]);
                    statement.execute();
                }
            }
        } catch (Exception $e) {
            //
        }
    }


}

*/
