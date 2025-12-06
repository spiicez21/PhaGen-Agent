"""
Test script for authentication endpoints.
Run this after starting the backend server to verify auth is working.
"""

import requests
import json
from datetime import datetime

BASE_URL = "http://localhost:8000/api"

def print_section(title):
    print("\n" + "="*60)
    print(f"  {title}")
    print("="*60)

def test_signup():
    """Test user signup"""
    print_section("TEST 1: User Signup")
    
    test_email = f"test_{datetime.now().timestamp()}@example.com"
    payload = {
        "email": test_email,
        "name": "Test User",
        "password": "testpass123"
    }
    
    print(f"\n📤 Sending signup request...")
    print(f"   Email: {test_email}")
    
    response = requests.post(f"{BASE_URL}/auth/signup", json=payload)
    
    if response.status_code == 201:
        print("✅ Signup successful!")
        user_data = response.json()
        print(f"   User ID: {user_data['id']}")
        print(f"   Email: {user_data['email']}")
        print(f"   Name: {user_data['name']}")
        return test_email, payload["password"]
    else:
        print(f"❌ Signup failed: {response.status_code}")
        print(f"   Error: {response.text}")
        return None, None

def test_login(email, password):
    """Test user login"""
    print_section("TEST 2: User Login")
    
    payload = {
        "email": email,
        "password": password
    }
    
    print(f"\n📤 Sending login request...")
    print(f"   Email: {email}")
    
    response = requests.post(f"{BASE_URL}/auth/login", json=payload)
    
    if response.status_code == 200:
        print("✅ Login successful!")
        token_data = response.json()
        token = token_data["access_token"]
        print(f"   Token type: {token_data['token_type']}")
        print(f"   Access token: {token[:50]}...")
        return token
    else:
        print(f"❌ Login failed: {response.status_code}")
        print(f"   Error: {response.text}")
        return None

def test_get_current_user(token):
    """Test getting current user info"""
    print_section("TEST 3: Get Current User")
    
    headers = {
        "Authorization": f"Bearer {token}"
    }
    
    print(f"\n📤 Sending authenticated request...")
    
    response = requests.get(f"{BASE_URL}/auth/me", headers=headers)
    
    if response.status_code == 200:
        print("✅ Retrieved user info successfully!")
        user_data = response.json()
        print(f"   User ID: {user_data['id']}")
        print(f"   Email: {user_data['email']}")
        print(f"   Name: {user_data['name']}")
        print(f"   Active: {user_data['is_active']}")
        return True
    else:
        print(f"❌ Failed to get user: {response.status_code}")
        print(f"   Error: {response.text}")
        return False

def test_invalid_login():
    """Test login with invalid credentials"""
    print_section("TEST 4: Invalid Login (Expected to Fail)")
    
    payload = {
        "email": "nonexistent@example.com",
        "password": "wrongpassword"
    }
    
    print(f"\n📤 Sending login request with invalid credentials...")
    
    response = requests.post(f"{BASE_URL}/auth/login", json=payload)
    
    if response.status_code == 401:
        print("✅ Correctly rejected invalid credentials!")
        return True
    else:
        print(f"❌ Unexpected response: {response.status_code}")
        return False

def test_unauthorized_access():
    """Test accessing protected endpoint without token"""
    print_section("TEST 5: Unauthorized Access (Expected to Fail)")
    
    print(f"\n📤 Sending request without token...")
    
    response = requests.get(f"{BASE_URL}/auth/me")
    
    if response.status_code == 403:
        print("✅ Correctly rejected unauthorized access!")
        return True
    else:
        print(f"❌ Unexpected response: {response.status_code}")
        return False

def test_create_job_authenticated(token):
    """Test creating a job with authentication"""
    print_section("TEST 6: Create Job (Authenticated)")
    
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json"
    }
    
    payload = {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "molecule": "Aspirin"
    }
    
    print(f"\n📤 Creating job with authentication...")
    print(f"   Molecule: {payload['molecule']}")
    
    response = requests.post(f"{BASE_URL}/jobs", headers=headers, json=payload)
    
    if response.status_code in [200, 201]:
        print("✅ Job created successfully!")
        job_data = response.json()
        print(f"   Job ID: {job_data.get('job_id', 'N/A')}")
        print(f"   Status: {job_data.get('status', 'N/A')}")
        return True
    else:
        print(f"⚠️  Job creation response: {response.status_code}")
        print(f"   Note: This may fail if job creation requires additional setup")
        return False

def main():
    print("\n🧪 PhaGen Authentication System Test Suite")
    print("=" * 60)
    print("Testing backend at:", BASE_URL)
    print("=" * 60)
    
    # Test 1: Signup
    email, password = test_signup()
    if not email:
        print("\n❌ Signup test failed. Cannot continue.")
        return
    
    # Test 2: Login
    token = test_login(email, password)
    if not token:
        print("\n❌ Login test failed. Cannot continue.")
        return
    
    # Test 3: Get current user
    test_get_current_user(token)
    
    # Test 4: Invalid login
    test_invalid_login()
    
    # Test 5: Unauthorized access
    test_unauthorized_access()
    
    # Test 6: Create authenticated job
    test_create_job_authenticated(token)
    
    # Summary
    print_section("TEST SUMMARY")
    print("\n✅ All authentication tests completed!")
    print("\nTest credentials:")
    print(f"   Email: {email}")
    print(f"   Password: {password}")
    print(f"   Token: {token[:50]}...")
    print("\nYou can use these credentials to test the frontend.")
    print("\n" + "="*60)

if __name__ == "__main__":
    try:
        main()
    except requests.exceptions.ConnectionError:
        print("\n❌ ERROR: Could not connect to backend server")
        print("   Make sure the backend is running on http://localhost:8000")
        print("\n   Start it with:")
        print("   cd backend")
        print("   python -m uvicorn app.main:app --reload")
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
