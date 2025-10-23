import React, { useState, useEffect } from 'react';
import Header from './components/Header';
import Hero from './components/Hero';
import Projects from './components/Projects';
import Footer from './components/Footer';
import ProjectPage from './pages/ProjectPage';

interface View {
  page: 'home' | 'project';
  projectId?: string | null;
}

const HomePage: React.FC<{ onNavigate: (page: 'project', projectId: string) => void }> = ({ onNavigate }) => (
  <>
    <Hero />
    <main className="container mx-auto px-4">
      <Projects onNavigate={onNavigate} />
    </main>
  </>
);

const App: React.FC = () => {
  const [view, setView] = useState<View>({ page: 'home', projectId: null });

  const handleNavigate = (page: 'home' | 'project', projectId: string | null = null) => {
    setView({ page, projectId });
    window.scrollTo(0, 0);
  };
  
  const handleScrollTo = (id: string) => {
    if (view.page !== 'home') {
      handleNavigate('home');
      // Wait for the home page to render before scrolling
      setTimeout(() => {
        document.getElementById(id)?.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }, 100);
    } else {
      document.getElementById(id)?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  };

  return (
    <div className="bg-green-800 text-green-200 min-h-screen font-sans">
      <Header onNavigateHome={() => handleNavigate('home')} onScrollTo={handleScrollTo} />
      {view.page === 'home' && <HomePage onNavigate={handleNavigate} />}
      {view.page === 'project' && view.projectId && (
        <ProjectPage projectId={view.projectId} onNavigateHome={() => handleNavigate('home')} onScrollToProjects={() => handleScrollTo('projects')} />
      )}
      <Footer />
    </div>
  );
};

export default App;
