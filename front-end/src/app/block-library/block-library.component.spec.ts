import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';
import { MatTooltip, MatTooltipModule } from '@angular/material/tooltip';

import { BlockId } from '../block.interface';
import { BlockLibraryComponent } from './block-library.component';
import { BlockService } from '../block.service';
import { MockBlockService } from '../mock-block.service';
import { MockOutputService } from '../mock-output.service';
import { OutputService } from '../output.service';

describe('BlockLibraryComponent', () => {
  let component: BlockLibraryComponent;
  let fixture: ComponentFixture<BlockLibraryComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [
        BlockLibraryComponent
      ],
      imports: [
        MatCardModule,
        MatSnackBarModule,
        MatTooltipModule
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: OutputService, useClass: MockOutputService }
      ]
    });
    fixture = TestBed.createComponent(BlockLibraryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Adding Blocks From Library', () => {
    it ('should call blockService.addBlock with corresponding block ID when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const blockIDs: BlockId[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      blockIDs.forEach(blockID => {
        const button = fixture.debugElement.query(By.css(`#${blockID}`));
        button.triggerEventHandler('click', {});
        fixture.detectChanges();
        expect(blockService.addBlock).toHaveBeenCalledWith(blockID);
      });
      expect(blockService.addBlock).toHaveBeenCalledTimes(blockIDs.length);
    });

    it ('should be available when blocks are not being executed', () => {
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('#loaddata')).nativeElement.disabled).toEqual(false);
    });

    it ('should become disabled while blocks are being executed', () => {
      component.executingBlocks = true;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const blockIDs: string[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).nativeElement.disabled).toEqual(true);
      });
    });

    it ('should become available once blocks have stopped being executed', () => {
      component.executingBlocks = true;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('#loaddata')).nativeElement.disabled).toEqual(false);
    });

    it ('should only allow Load Data block to be added when canvas is empty', () => {
      component.blockList = [];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const availableBlocks: string[] = ['loaddata'];
      const unavailableBlocks: string[] = ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      availableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).nativeElement.disabled).toEqual(false);
      });
      unavailableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).nativeElement.disabled).toEqual(true);
      });
    });

    it ('should only allow possible child blocks of the last block on the canvas to be added when canvas is not empty', () => {
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'],
        parameters: [],
      }];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const availableBlocks: string[] = ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'];
      const unavailableBlocks: string[] = ['loaddata', 'pca', 'integration', 'runumap'];
      availableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).nativeElement.disabled).toEqual(false);
      });
      unavailableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).nativeElement.disabled).toEqual(true);
      });
    });

    it ('should have a tooltip enabled on the add button only when a block cannot be added', () => {
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'],
        parameters: [],
      }];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const availableBlocks: string[] = ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'];
      const unavailableBlocks: string[] = ['loaddata', 'pca', 'integration', 'runumap'];
      availableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).disabled).toEqual(true);
      });
      unavailableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).disabled).toEqual(false);
      });
    });

    it ('should have an appropriate tooltip message when blocks are being executed', () => {
      component.blockList = [];
      component.executingBlocks = true;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const blockIDs: BlockId[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).disabled).toEqual(false);
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).message).toEqual('Blocks cannot be added while blocks are being executed.');
      });
    });

    it ('should have an appropriate tooltip message when no blocks have been added', () => {
      component.blockList = [];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const unavailableBlocks: string[] = ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      unavailableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).message).toEqual('This block cannot be added to an empty canvas.');
      });
    });

    it ('should have an appropriate tooltip message when a block cannot follow the existing final block (which starts with a vowel)', () => {
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'A Vowel',
        possibleChildBlocks: ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'],
        parameters: [],
      }];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const unavailableBlocks: string[] = ['loaddata', 'pca', 'integration', 'runumap'];
      unavailableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).message).toEqual('This block cannot be immediately below an A Vowel block.');
      });
    });

    it ('should have an appropriate tooltip message when a block cannot follow the existing final block (which does not start with a vowel)', () => {
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Not A Vowel',
        possibleChildBlocks: ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'],
        parameters: [],
      }];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const unavailableBlocks: string[] = ['loaddata', 'pca', 'integration', 'runumap'];
      unavailableBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css(`#${blockID}`)).injector.get<MatTooltip>(MatTooltip).message).toEqual('This block cannot be immediately below a Not A Vowel block.');
      });
    });
  });

  describe('Help text for blocks', () => {
    it ('should have an appropriate tooltip message for each block when user hovers over icon', () => {
      interface BlockHelpText {
        id: BlockId
        text: string
      }
      const blockHelpTexts: BlockHelpText[] = [
        {id: 'loaddata', text: 'This step allows us to read our data into SCAMPI ready for analysis!'},
        {id: 'basicfiltering', text: 'Remove cells and genes that are not expressing much'},
        {id: 'qcplots', text: 'Plot some metrics to see if we have any outliers in our data'},
        {id: 'qcfiltering', text: 'Apply a threshold to remove outliers that could be doublets or stressed/dying cells'},
        {id: 'variablegenes', text: 'This selects a set of genes that explain the most variation in our data and explain most of the underlying biology'},
        {id: 'pca', text: 'PCA helps us group genes together in principle components that explain the biological differences in our data'},
        {id: 'integration', text: 'Unwanted sources of variation need to be corrected for, revealing differences that are driven by the biology and not batch effects'},
        {id: 'runumap', text: 'This allows us to visualise our complex data in a simpler way'}
      ];
      blockHelpTexts.forEach(block => {
        expect(fixture.debugElement.query(By.css(`#${block.id}-help`)).injector.get<MatTooltip>(MatTooltip).disabled).toEqual(false);
        expect(fixture.debugElement.query(By.css(`#${block.id}-help`)).injector.get<MatTooltip>(MatTooltip).message).toEqual(block.text);
      });
    });
  });
});
