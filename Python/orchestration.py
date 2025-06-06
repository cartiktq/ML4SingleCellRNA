if __name__ == "__main__":
    steps = [
        LoadData(output_path="step1"),
        PreprocessingQC(input_path="step1", output_path="step2"),
        DimensionalityReduction(input_path="step2", output_path="step3"),
        ClusteringStratification(input_path="step3", output_path="step4"),
        BiomarkerIdentification(input_path="step4", output_path="step5"),
        TrajectoryAnalysis(input_path="step5", output_path="step6"),
        Visualization(input_path="step6", output_path="step7")
    ]

    for step in steps:
        print(f"Running: {step.__class__.__name__}")
        step.run()
