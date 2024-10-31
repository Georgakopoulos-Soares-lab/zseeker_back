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
	log.Printf("Executing command: %s", cmd.String())
	output, err := cmd.CombinedOutput()
	if err != nil {
		log.Printf("Error executing ZSeeker: %v", err)
		log.Printf("Command output: %s", string(output))
		log.Printf("Command error details: %v", err)
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to execute ZSeeker", "details": string(output)})
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
	params["GT_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("GT_weight", "1.0"), 64)
	params["AC_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("AC_weight", "1.0"), 64)
	params["AT_weight"], _ = strconv.ParseFloat(c.DefaultPostForm("AT_weight", "0.5"), 64)
	params["mismatch_penalty_starting_value"], _ = strconv.Atoi(c.DefaultPostForm("mismatch_penalty_starting_value", "1"))
	params["mismatch_penalty_linear_delta"], _ = strconv.Atoi(c.DefaultPostForm("mismatch_penalty_linear_delta", "2"))
	params["mismatch_penalty_type"] = c.DefaultPostForm("mismatch_penalty_type", "linear")
	
	// Parse consecutive_AT_scoring as a slice of floats
	consecutiveATScoringStr := c.DefaultPostForm("consecutive_AT_scoring", "0.5,0.5,0.5,0.5,0.0,0.0,-5.0,-100.0")
	consecutiveATScoringStrs := strings.Split(consecutiveATScoringStr, ",")
	var consecutiveATScoring []float64
	for _, v := range consecutiveATScoringStrs {
		score, err := strconv.ParseFloat(strings.TrimSpace(v), 64)
		if err == nil {
			consecutiveATScoring = append(consecutiveATScoring, score)
		} else {
			// Handle conversion error (optional)
			fmt.Printf("Warning: Could not parse '%s' as float: %v\n", v, err)
		}
	}
	
	params["consecutive_AT_scoring"] = consecutiveATScoring  // Store the slice directly
	params["n_jobs"], _ = strconv.Atoi(c.DefaultPostForm("n_jobs", "8"))
	params["cadence_reward"], _ = strconv.ParseFloat(c.DefaultPostForm("cadence_reward", "1.0"), 64)
	params["method"] = c.DefaultPostForm("method", "transitions")
	params["threshold"], _ = strconv.Atoi(c.DefaultPostForm("threshold", "50"))

	return params
}

func prepareCommand(params map[string]interface{}, filePath string) *exec.Cmd {
	cmdArgs := []string{
		"--path", filePath,
		"--GC_weight", fmt.Sprintf("%v", params["GC_weight"]),
		"--AT_weight", fmt.Sprintf("%v", params["AT_weight"]),
		"--GT_weight", fmt.Sprintf("%v", params["GT_weight"]),
		"--AC_weight", fmt.Sprintf("%v", params["AC_weight"]),
		"--mismatch_penalty_starting_value", fmt.Sprintf("%v", params["mismatch_penalty_starting_value"]),
		"--mismatch_penalty_linear_delta", fmt.Sprintf("%v", params["mismatch_penalty_linear_delta"]),
		"--mismatch_penalty_type", params["mismatch_penalty_type"].(string),
		"--method", params["method"].(string),
		"--n_jobs", fmt.Sprintf("%v", params["n_jobs"]),
		"--threshold", fmt.Sprintf("%v", params["threshold"]),
	}

	// Handle consecutive_AT_scoring as a list of floats
	consecutiveATScoring := params["consecutive_AT_scoring"].([]float64)
	for _, score := range consecutiveATScoring {
		cmdArgs = append(cmdArgs, "--consecutive_AT_scoring", fmt.Sprintf("%v", score))
	}

	return exec.Command("ZSeeker", cmdArgs...)
}

func saveFile(fileHeader *multipart.FileHeader, destination string) error {
	log.Println("Saving the uploaded file to disk...")
	return nil // Placeholder implementation for saving file logic
}
