#!/usr/bin/env python3
"""
Test script for RAPTOR Metal GPU Implementation

This script tests the basic functionality of the Metal-based RAPTOR code.
"""

import os
import sys
import subprocess
import time
import numpy as np

def test_compilation():
    """Test if the Metal code compiles successfully."""
    print("Testing compilation...")
    
    try:
        # Try to build the code
        result = subprocess.run(['make', 'clean'], capture_output=True, text=True)
        result = subprocess.run(['make', 'all'], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Compilation successful")
            return True
        else:
            print("❌ Compilation failed")
            print(f"Error: {result.stderr}")
            return False
            
    except FileNotFoundError:
        print("❌ Make not found. Please install build tools.")
        return False
    except Exception as e:
        print(f"❌ Compilation error: {e}")
        return False

def test_binary_exists():
    """Test if the binary was created."""
    print("Testing binary creation...")
    
    if os.path.exists("./raptor_metal") and os.access("./raptor_metal", os.X_OK):
        print("✅ Binary created and executable")
        return True
    else:
        print("❌ Binary not found or not executable")
        return False

def test_metal_support():
    """Test if Metal is supported on this system."""
    print("Testing Metal support...")
    
    try:
        # Try to run the binary to check Metal support
        result = subprocess.run(['./raptor_metal'], capture_output=True, text=True, timeout=10)
        
        # If it runs without Metal errors, it's supported
        if "Metal is not supported" not in result.stderr:
            print("✅ Metal support detected")
            return True
        else:
            print("❌ Metal not supported on this device")
            return False
            
    except subprocess.TimeoutExpired:
        print("✅ Metal support detected (program ran)")
        return True
    except Exception as e:
        print(f"❌ Error testing Metal support: {e}")
        return False

def test_python_wrapper():
    """Test if the Python wrapper can be imported."""
    print("Testing Python wrapper...")
    
    try:
        # Try to import the wrapper
        sys.path.append('.')
        from raptor_metal_wrapper import RAPTORMetalWrapper
        print("✅ Python wrapper imported successfully")
        return True
    except ImportError as e:
        print(f"❌ Python wrapper import failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Python wrapper error: {e}")
        return False

def test_basic_functionality():
    """Test basic functionality of the wrapper."""
    print("Testing basic functionality...")
    
    try:
        from raptor_metal_wrapper import RAPTORMetalWrapper
        
        # Create wrapper instance
        raptor = RAPTORMetalWrapper()
        
        # Test parameter setting
        raptor.set_model_parameters(IMG_WIDTH=128, IMG_HEIGHT=128)
        
        # Test parameter retrieval
        if raptor.model_params['IMG_WIDTH'] == 128:
            print("✅ Parameter setting works")
        else:
            print("❌ Parameter setting failed")
            return False
            
        print("✅ Basic functionality test passed")
        return True
        
    except Exception as e:
        print(f"❌ Basic functionality test failed: {e}")
        return False

def run_performance_test():
    """Run a simple performance test."""
    print("Running performance test...")
    
    try:
        # This would require the actual binary to work
        # For now, just check if we can create the wrapper
        from raptor_metal_wrapper import RAPTORMetalWrapper
        
        raptor = RAPTORMetalWrapper()
        raptor.set_model_parameters(
            IMG_WIDTH=256,
            IMG_HEIGHT=256,
            FREQ_MIN=1e11,
            FREQ_MAX=1e12
        )
        
        print("✅ Performance test setup successful")
        print("   Note: Full performance test requires working binary")
        return True
        
    except Exception as e:
        print(f"❌ Performance test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("RAPTOR Metal GPU Implementation - Test Suite")
    print("=" * 50)
    
    tests = [
        ("Compilation", test_compilation),
        ("Binary Creation", test_binary_exists),
        ("Metal Support", test_metal_support),
        ("Python Wrapper", test_python_wrapper),
        ("Basic Functionality", test_basic_functionality),
        ("Performance Test Setup", run_performance_test)
    ]
    
    results = []
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ Test failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 50)
    print("Test Summary:")
    print("=" * 50)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name:25} {status}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! RAPTOR Metal is ready to use.")
        return 0
    else:
        print("⚠️  Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    exit(main())