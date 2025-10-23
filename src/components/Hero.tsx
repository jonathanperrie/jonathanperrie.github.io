import React from 'react';
import { SOCIAL_LINKS } from '../constants';
import GitHubIcon from './icons/GitHubIcon';
import LinkedInIcon from './icons/LinkedInIcon';
import GoogleScholarIcon from './icons/GoogleScholarIcon';
import DocumentIcon from './icons/DocumentIcon';

const Hero: React.FC = () => {
  return (
    <section id="home" className="h-[calc(100vh-80px)] min-h-[600px] flex items-center justify-center text-center">
      <div className="container mx-auto px-4">
        <h1 className="text-4xl md:text-6xl font-extrabold text-green-100 mb-4">
          Hi, I'm a Senior Frontend Engineer
        </h1>
        <p className="text-xl md:text-2xl text-green-400 mb-8 max-w-3xl mx-auto">
          I build beautiful, responsive, and high-performance web applications with a focus on user experience.
        </p>
        <div className="flex justify-center items-center space-x-6">
          <a href={SOCIAL_LINKS.linkedIn} aria-label="LinkedIn" className="text-green-400 hover:text-lime-400 transition-colors">
            <LinkedInIcon className="w-8 h-8" />
          </a>
          <a href={SOCIAL_LINKS.github} aria-label="GitHub" className="text-green-400 hover:text-lime-400 transition-colors">
            <GitHubIcon className="w-8 h-8" />
          </a>
          <a href={SOCIAL_LINKS.googleScholar} aria-label="Google Scholar" className="text-green-400 hover:text-lime-400 transition-colors">
            <GoogleScholarIcon className="w-8 h-8" />
          </a>
          <a href={SOCIAL_LINKS.resume} aria-label="Resume" className="text-green-400 hover:text-lime-400 transition-colors">
            <DocumentIcon className="w-8 h-8" />
          </a>
        </div>
      </div>
    </section>
  );
};

export default Hero;
