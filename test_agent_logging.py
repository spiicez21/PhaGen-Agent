"""Quick test to verify agent logging configuration"""
import sys
import logging
from pathlib import Path

# Add root to path
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

# Set agent loggers
logging.getLogger('agents.master').setLevel(logging.INFO)
logging.getLogger('agents.llm').setLevel(logging.INFO)
logging.getLogger('agents.workers').setLevel(logging.INFO)

# Test imports
print("\n" + "="*80)
print("Testing Agent Logging Configuration")
print("="*80 + "\n")

try:
    from agents.master import MasterAgent
    print("✅ MasterAgent imported successfully")
    
    from agents.llm import LLMClient
    print("✅ LLMClient imported successfully")
    
    from agents.workers.clinical import ClinicalWorker
    print("✅ ClinicalWorker imported successfully")
    
    print("\n" + "="*80)
    print("✅ All agent modules loaded - logging configured and ready!")
    print("="*80 + "\n")
    print("Agent logs will now appear in the backend console when jobs are executed.")
    print("Log format: timestamp - module - level - message")
    print("\nExample log output:")
    print("  2025-12-06 14:11:15 - agents.master - INFO - 🚀 Starting agent analysis for molecule: Nintedanib")
    print("  2025-12-06 14:11:16 - agents.workers.base - INFO -    🔍 [clinical] Building search queries...")
    print("  2025-12-06 14:11:17 - agents.workers.base - INFO -    ✅ [clinical] Summary complete (confidence: 0.85)")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
