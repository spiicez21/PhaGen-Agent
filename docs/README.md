# PhaGen Documentation

Welcome to the PhaGen documentation. This directory contains comprehensive documentation for the PhaGen Agentic system, a pharmaceutical evidence aggregation and innovation assessment platform.

## 📚 Documentation Structure

The documentation is organized in a logical sequence to help you understand, build, and deploy PhaGen:

### Core Documentation

1. **[01-ARCHITECTURE.md](./01-ARCHITECTURE.md)**
   - System architecture overview
   - Component interactions and data flow
   - Technology stack and design decisions

2. **[02-GETTING-STARTED.md](./02-GETTING-STARTED.md)**
   - Setup and installation instructions
   - Development environment configuration
   - Quick start guide for running PhaGen locally

3. **[03-PIPELINE-OVERVIEW.md](./03-PIPELINE-OVERVIEW.md)**
   - Data processing pipeline
   - Evidence retrieval and aggregation workflow
   - Worker agent execution flow

4. **[04-UI-WIREFRAMES.md](./04-UI-WIREFRAMES.md)**
   - User interface design specifications
   - Component layouts and user flows
   - Frontend implementation guidelines

5. **[05-PROJECT-ROADMAP.md](./05-PROJECT-ROADMAP.md)**
   - Project execution tracker
   - Completed features and milestones
   - Pending tasks and future enhancements
   - Implementation priorities

### Feature Documentation

6. **[06-ML-FEATURES.md](./06-ML-FEATURES.md)**
   - Machine learning capabilities
   - Disease mapping model
   - Reranker fine-tuning pipeline
   - Custom embedding training
   - Active learning implementation

7. **[07-ADVANCED-FEATURES.md](./07-ADVANCED-FEATURES.md)**
   - Patent claim-level semantic matching
   - Continuous crawler with change detection
   - Paid market API integration
   - Multi-region regulatory rule engine
   - Advanced analytics and workflows

8. **[08-ENTERPRISE-SECURITY.md](./08-ENTERPRISE-SECURITY.md)**
   - Hardware-backed key storage (HSM)
   - Signed container images and supply-chain verification
   - Air-gapped deployment mode
   - Security architecture and compliance
   - Production deployment checklists

9. **[09-IMPLEMENTATION-STATUS.md](./09-IMPLEMENTATION-STATUS.md)**
   - Current implementation status
   - Completed phases and features
   - Testing and validation status
   - Production readiness assessment

## 🎯 Quick Navigation

### For New Developers
Start with:
1. [02-GETTING-STARTED.md](./02-GETTING-STARTED.md) - Setup your environment
2. [01-ARCHITECTURE.md](./01-ARCHITECTURE.md) - Understand the system
3. [03-PIPELINE-OVERVIEW.md](./03-PIPELINE-OVERVIEW.md) - Learn the data flow

### For Product Managers
Focus on:
1. [05-PROJECT-ROADMAP.md](./05-PROJECT-ROADMAP.md) - Track progress and priorities
2. [04-UI-WIREFRAMES.md](./04-UI-WIREFRAMES.md) - Review UI/UX specifications
3. [09-IMPLEMENTATION-STATUS.md](./09-IMPLEMENTATION-STATUS.md) - Check feature completeness

### For DevOps/Security
Review:
1. [08-ENTERPRISE-SECURITY.md](./08-ENTERPRISE-SECURITY.md) - Security infrastructure
2. [02-GETTING-STARTED.md](./02-GETTING-STARTED.md) - Deployment procedures
3. [01-ARCHITECTURE.md](./01-ARCHITECTURE.md) - Infrastructure requirements

### For Data Scientists
Explore:
1. [06-ML-FEATURES.md](./06-ML-FEATURES.md) - ML models and training
2. [07-ADVANCED-FEATURES.md](./07-ADVANCED-FEATURES.md) - Advanced analytics
3. [03-PIPELINE-OVERVIEW.md](./03-PIPELINE-OVERVIEW.md) - Data processing

## 📖 Additional Resources

### External Documentation
- **Infra Scripts**: See `../infra/scripts/SIGNING_GUIDE.md` for container signing
- **Air-Gap Deployment**: See `../infra/airgap/AIRGAP_DEPLOYMENT_GUIDE.md` for offline deployment
- **API Documentation**: Run the backend and visit `/docs` for FastAPI Swagger UI
- **Frontend Components**: See `../frontend/README.md` for React component docs

### Code Examples
- **Backend**: `../backend/` - FastAPI application and worker agents
- **Frontend**: `../frontend/` - Next.js dashboard application
- **Crawler**: `../crawler/` - Crawlee data acquisition scripts
- **Indexes**: `../indexes/` - Vector database and embedding pipelines
- **Infrastructure**: `../infra/` - Docker Compose, Kubernetes, security tools

## 🔄 Documentation Updates

This documentation is actively maintained. When making changes to the system:

1. **Update relevant documentation files** to reflect new features or changes
2. **Follow the numbering convention** for new major documentation sections
3. **Update this README.md** if you add new documentation files
4. **Keep code examples current** with actual implementation
5. **Document breaking changes** prominently in relevant sections

## 📋 Version Information

- **Last Updated**: December 2024
- **PhaGen Version**: 1.0.0
- **Documentation Version**: 1.0.0

## 🤝 Contributing

When contributing to documentation:

- Use clear, concise language
- Include code examples where applicable
- Add diagrams for complex workflows (use Mermaid syntax)
- Cross-reference related documentation sections
- Keep technical depth appropriate for target audience
- Test all commands and code snippets before documenting

## 📧 Support

For questions about this documentation:
- Open an issue in the repository
- Contact the development team
- Review the implementation files for code-level details

---

**Happy Building!** 🚀
