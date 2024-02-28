import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockService } from '../block.service';
import { CodeBlockComponent } from './code-block.component';

describe('CodeBlockComponent', () => {
  let component: CodeBlockComponent;
  let fixture: ComponentFixture<CodeBlockComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CodeBlockComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
      ],
      providers: [BlockService],
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Removing Blocks', () => {
    it ('blockService.removeBlock should be called with loaddata when remove button is clicked', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('loaddata');
    });

    it ('blockService.removeBlock should be called with basicfiltering when remove button is clicked', () => {
      component.block = {
        blockId: 'basicfiltering',
        title: 'Basic Filtering',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('basicfiltering');
    });

    it ('blockService.removeBlock should be called with qcplots when remove button is clicked', () => {
      component.block = {
        blockId: 'qcplots',
        title: 'Quality Control Plots',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('qcplots');
    });

    it ('blockService.removeBlock should be called with qcfiltering when remove button is clicked', () => {
      component.block = {
        blockId: 'qcfiltering',
        title: 'Quality Control Filtering',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('qcfiltering');
    });

    it ('blockService.removeBlock should be called with variablegenes when remove button is clicked', () => {
      component.block = {
        blockId: 'variablegenes',
        title: 'Identify Highly Variable Genes',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('variablegenes');
    });

    it ('blockService.removeBlock should be called with pca when remove button is clicked', () => {
      component.block = {
        blockId: 'pca',
        title: 'Principle Component Analysis',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('pca');
    });
  });
});
