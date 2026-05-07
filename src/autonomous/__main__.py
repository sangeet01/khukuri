import argparse
from .controller import AutonomousLoopController, Mode

def main():
    parser = argparse.ArgumentParser(description="Khukuri Autonomous Loop CLI")
    parser.add_argument("--mode", type=str, choices=["supervised", "lights_out"], default="supervised")
    parser.add_argument("--project", type=str, default="default")
    parser.add_argument("--pathogen", type=str, default="S. aureus MRSA")
    
    args = parser.parse_args()
    
    mode = Mode.LIGHTS_OUT if args.mode == "lights_out" else Mode.SUPERVISED
    
    controller = AutonomousLoopController(
        mode=mode,
        project=args.project,
        pathogen=args.pathogen
    )
    
    controller.run_cycle()

if __name__ == "__main__":
    main()
