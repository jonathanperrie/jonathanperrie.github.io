import React from 'react';

interface HeaderProps {
  onNavigateHome: () => void;
  onScrollTo: (id: string) => void;
}

const Header: React.FC<HeaderProps> = ({ onNavigateHome, onScrollTo }) => {
  const handleNavClick = (e: React.MouseEvent<HTMLButtonElement>, id: string) => {
    e.preventDefault();
    onScrollTo(id);
  };

  return (
    <header className="sticky top-0 bg-green-800/80 backdrop-blur-sm z-50">
      <nav className="container mx-auto px-4 py-4 flex justify-between items-center">
        <button onClick={onNavigateHome} className="text-2xl font-bold text-lime-400 hover:text-lime-300 transition-colors">
          Hello
        </button>
        <ul className="flex space-x-6">
          <li>
            <button onClick={(e) => handleNavClick(e, 'projects')} className="text-lg text-green-300 hover:text-lime-400 transition-colors">
              Projects
            </button>
          </li>
        </ul>
      </nav>
      <div className="h-px bg-gradient-to-r from-transparent via-lime-500 to-transparent"></div>
    </header>
  );
};

export default Header;