class WorkflowStep:
    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path

    def run(self):
        raise NotImplementedError("Each subclass must implement the run() method.")
