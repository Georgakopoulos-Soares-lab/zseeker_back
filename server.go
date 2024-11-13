package main

import (
	"fmt"
	"log"
	"mime/multipart"
	"net/http"
	"os"
	"os/exec"
	"encoding/csv"
	"path/filepath"
	"strconv"
	"strings"
	"github.com/gin-gonic/gin"
)

func main() {
	
	// Initialize a new Gin router with default middleware: logger and recovery middleware
	r := gin.Default()
	
	// Setup CORS Middleware
	r.Use(func(c *gin.Context) {
		c.Writer.Header().Set("Access-Control-Allow-Origin", "*") // Allow all origins
		c.Writer.Header().Set("Access-Control-Allow-Methods", "GET, POST, PUT, DELETE, OPTIONS")
		c.Writer.Header().Set("Access-Control-Allow-Headers", "Content-Type, Authorization")
		if c.Request.Method == "OPTIONS" {
			c.AbortWithStatus(http.StatusNoContent)
			return
		}
		c.Next()
	})
	
	// Routes
	r.POST("/submit-job", submitJob)
	
	// Start Server
	log.Println("Starting the server on port 8080...")
	r.Run(":8080")
}

func submitJob(c *gin.Context) {
	log.Println("Handling /submit-job request...")

	file, err := c.FormFile("fasta")
	if err != nil {
		log.Printf("Failed to read the uploaded file: %v", err)
		c.JSON(http.StatusBadRequest, gin.H{"error": "Failed to read the fasta file"})
		return
	}

	log.Println("Saving uploaded file...")
	filePath := filepath.Join("./uploads", filepath.Base(file.Filename))
	err = c.SaveUploadedFile(file, filePath)
	if err != nil {
		log.Printf("Failed to save file: %v", err)
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to save the file"})
		return
	}

	params := parseParams(c)
	if params == nil {
		c.JSON(http.StatusBadRequest, gin.H{"error": "Invalid input parameters"})
		return
	}

	log.Println("Preparing to execute ZSeeker...")
	cmd := prepareCommand(params, filePath)
	log.Printf("Executing command: ZSeeker %s", strings.Join(cmd.Args[1:], " "))
	output, err := cmd.CombinedOutput()
	log.Printf("Command output:\n%s", string(output))
	if err != nil {
		log.Printf("Error executing ZSeeker: %v", err)
		log.Printf("Command error details: %v", err)
		c.JSON(http.StatusInternalServerError, gin.H{
			"error": "Failed to execute ZSeeker",
			"details": string(output),
			"command": fmt.Sprintf("ZSeeker %s", strings.Join(cmd.Args[1:], " ")),
		})
		return
	}

	log.Println("Command executed successfully. Reading Results")

	dat, err := os.ReadFile("./extractions_zdna_human/input_zdna_score.csv")
    if err != nil {
		log.Printf("Failed to read input file %v", err)
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to read results file"})
		return
	} 
    // Parse the CSV file
	csvReader := csv.NewReader(strings.NewReader(string(dat)))
	records, err := csvReader.ReadAll()
	if err != nil {
		log.Printf("Failed to parse CSV file %v", err)
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to parse CSV file"})
		return
	}

	// Separate headers and data
	csvHeaders := records[0]
	csvData := records[1:]

	// Log and respond with headers and data in JSON format
	log.Println("CSV Headers:", csvHeaders)
	log.Println("CSV Data:", csvData)
	c.JSON(http.StatusOK, gin.H{
		"message":     "Job completed successfully",
		"csv_headers": csvHeaders,
		"csv_data":    csvData,
	})
}

func parseParams(c *gin.Context) map[string]interface{} {
    params := make(map[string]interface{})

    params["GC_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("GC_weight", "7.0"), 64)
    params["GT_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("GT_weight", "1.25"), 64)
    params["AC_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("AC_weight", "1.25"), 64)
    params["AT_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("AT_weight", "0.5"), 64)
    params["mismatch_penalty_starting_value"], _ = strconv.Atoi(c.DefaultPostForm("mismatch_penalty_starting_value", "3"))
    params["mismatch_penalty_linear_delta"], _ = strconv.Atoi(c.DefaultPostForm("mismatch_penalty_linear_delta", "3"))
    params["mismatch_penalty_type"] = c.DefaultPostForm("mismatch_penalty_type", "linear")
    
    // Parse consecutive_AT_scoring as a slice of floats
    defaultScoring := "0.5,0.5,0.5,0.5,0.0,0.0,-5.0,-100.0"
    consecutiveATScoringStr := strings.TrimSpace(c.DefaultPostForm("consecutive_AT_scoring", defaultScoring))
    consecutiveATScoringStrs := strings.Split(consecutiveATScoringStr, ",")
    consecutiveATScoring := make([]float64, 0, len(consecutiveATScoringStrs))
    
    for _, v := range consecutiveATScoringStrs {
        score, err := strconv.ParseFloat(strings.TrimSpace(v), 64)
        if err == nil {
            consecutiveATScoring = append(consecutiveATScoring, score)
        } else {
            log.Printf("Warning: Could not parse '%s' as float: %v\n", v, err)
        }
    }
    
    // If parsing failed completely, use default values
    if len(consecutiveATScoring) == 0 {
        defaultValues := []float64{0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0}
        consecutiveATScoring = defaultValues
    }
    
    params["consecutive_AT_scoring"] = consecutiveATScoring
    params["n_jobs"], _ = strconv.Atoi(c.DefaultPostForm("n_jobs", "8"))
    params["threshold"], _ = strconv.Atoi(c.DefaultPostForm("threshold", "50"))

    return params
}


func prepareCommand(params map[string]interface{}, filePath string) *exec.Cmd {
    // Convert consecutive_AT_scoring array to the required format
    consecutiveATScoring := params["consecutive_AT_scoring"].([]float64)
    atScoringStr := strings.Builder{}
    for i, score := range consecutiveATScoring {
        if i > 0 {
            atScoringStr.WriteString(",")
        }
        atScoringStr.WriteString(fmt.Sprintf("%.1f", score)) // Keep %.1f for consecutive_AT_scoring
    }

    cmdArgs := []string{
        "--fasta", filePath,
        "--GC_weight", fmt.Sprintf("%.2f", params["GC_weight"]), // Changed to %.2f
        "--AT_weight", fmt.Sprintf("%.2f", params["AT_weight"]), // Changed to %.2f
        "--GT_weight", fmt.Sprintf("%.2f", params["GT_weight"]), // Changed to %.2f
        "--AC_weight", fmt.Sprintf("%.2f", params["AC_weight"]), // Changed to %.2f
        "--mismatch_penalty_starting_value", fmt.Sprintf("%d", params["mismatch_penalty_starting_value"]),
        "--mismatch_penalty_linear_delta", fmt.Sprintf("%d", params["mismatch_penalty_linear_delta"]),
        "--mismatch_penalty_type", params["mismatch_penalty_type"].(string),
        "--n_jobs", fmt.Sprintf("%d", params["n_jobs"]),
        "--threshold", fmt.Sprintf("%d", params["threshold"]),
        "--consecutive_AT_scoring", atScoringStr.String(),
        "--output_dir", "extractions_zdna_human",
    }

    cmd := exec.Command("ZSeeker", cmdArgs...)
    
    // Log the complete command for debugging
    log.Printf("Complete ZSeeker command: ZSeeker %s", strings.Join(cmdArgs, " "))
    
    return cmd
}


func saveFile(fileHeader *multipart.FileHeader, destination string) error {
	log.Println("Saving the uploaded file to disk...")
	return nil // Placeholder implementation for saving file logic
}
