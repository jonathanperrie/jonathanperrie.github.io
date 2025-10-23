import React, { useEffect, useState } from 'react';
import { PROJECTS } from '../constants';
import { Project } from '../types';

interface ProjectPageProps {
  projectId: string;
  onNavigateHome: () => void;
  onScrollToProjects: () => void;
}

const ProjectPage: React.FC<ProjectPageProps> = ({ projectId, onScrollToProjects }) => {
  const [project, setProject] = useState<Project | undefined>(undefined);

  useEffect(() => {
    const foundProject = PROJECTS.find(p => p.id === projectId);
    setProject(foundProject);
  }, [projectId]);

  if (!project) {
    return (
      <div className="container mx-auto px-4 py-20 text-center">
        <h2 className="text-2xl font-bold text-green-200">Project not found</h2>
        <button onClick={onScrollToProjects} className="text-lime-400 hover:text-lime-300 mt-4 inline-block">
          Back to all projects
        </button>
      </div>
    );
  }

  return (
    <main className="container mx-auto px-4 py-16 md:py-24">
      <div className="max-w-4xl mx-auto">
        <button onClick={onScrollToProjects} className="text-lime-400 hover:text-lime-300 mb-8 inline-block">
          Back to all projects
        </button>
        <h1 className="text-4xl md:text-5xl font-bold text-green-100 mb-4">{project.title}</h1>
        <div className="mb-8">
          {project.tags.map((tag, index) => (
            <span key={index} className="inline-block bg-lime-900/50 text-lime-300 text-xs font-semibold mr-2 mb-2 px-2.5 py-0.5 rounded-full">
              {tag}
            </span>
          ))}
        </div>
        <img src={project.imageUrl} alt={project.title} className="w-full h-auto rounded-lg shadow-lg mb-8" />
        <div className="text-lg text-green-300 space-y-6">
          <p>{project.detailedDescription}</p>
        </div>
        <div className="mt-12">
          <a 
            href={project.projectUrl} 
            className="inline-block bg-lime-500 text-green-900 font-bold py-3 px-6 rounded-lg hover:bg-lime-600 transition-colors"
          >
            View Project on GitHub
          </a>
        </div>
      </div>
    </main>
  );
};

export default ProjectPage;
